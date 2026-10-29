"""GPU integration of Chua's ODE and its variational equation using PyOpenCL."""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pyopencl as cl

CLASS_NAMES = {0: "escaped", 1: "equilibrium", 2: "nonpositive_or_near_zero",
               3: "uncertain", 4: "chaotic_candidate", 5: "numerical_failure"}
DEFAULT_SEEDS = [[0.7, 0.0, 0.0], [-0.2, 0.1, 0.0], [1.5, 0.02, -1.5]]


def rhs(t, state, parameters):
    """CPU reference vector field; parameters are alpha, beta, m0, m1."""
    a, b, m0, m1 = parameters
    x, y, z = state[:3]
    h = m1*x + 0.5*(m0-m1)*(abs(x+1)-abs(x-1))
    flow = [a*(y-x-h), x-y+z, -b*y]
    if len(state) == 6:
        u, v, w = state[3:]
        slope = m0 if abs(x)<1 else m1
        flow += [a*(v-(1+slope)*u), u-v+w, -b*v]
    return np.asarray(flow)


class GPU:
    def __init__(self, device_filter="", precision="float64"):
        devices = [d for p in cl.get_platforms() for d in p.get_devices()
                   if d.type & cl.device_type.GPU and device_filter.lower() in d.name.lower()]
        if not devices:
            raise RuntimeError("No matching OpenCL GPU. Install a GPU driver with OpenCL support. "
                               "This program does not silently fall back to CPU.")
        self.device = devices[0]
        if precision not in ("float32", "float64"):
            raise ValueError("precision must be float32 or float64")
        if precision == "float64" and not self.device.double_fp_config:
            raise RuntimeError("This GPU does not support float64; explicitly select float32.")
        self.dtype = np.dtype(precision)
        self.scalar = self.dtype.type
        self.context = cl.Context([self.device])
        self.queue = cl.CommandQueue(self.context, properties=cl.command_queue_properties.PROFILING_ENABLE)
        source = Path(__file__).with_name("kernels.cl").read_text(encoding="utf-8")
        self.program = cl.Program(self.context, source).build(
            options=["-DUSE_DOUBLE"] if precision == "float64" else [])
        self.advance = cl.Kernel(self.program, "advance")
        self.trace_kernel = cl.Kernel(self.program, "trace")
        self.info = {"name": self.device.name.strip(), "vendor": self.device.vendor,
                     "driver": self.device.driver_version, "opencl": self.device.version,
                     "precision": precision, "global_memory_bytes": self.device.global_mem_size,
                     "backend": "PyOpenCL", "fast_math": False}

    def buffer(self, array, readonly=False):
        flags = cl.mem_flags.READ_ONLY if readonly else cl.mem_flags.READ_WRITE
        return cl.Buffer(self.context, flags | cl.mem_flags.COPY_HOST_PTR, hostbuf=array)

    def scan(self, parameters, seeds=DEFAULT_SEEDS, dt=0.01, transient=300.0,
             duration=800.0, threshold=0.01, escape=100.0, chunk_steps=200,
             progress=None, initial_states=None):
        """Four independent time-window LLE estimates for each parameter/seed pair.

        A candidate must be bounded for the observation, have mean LLE > threshold,
        and have LLE > threshold/2 in *every* quarter. This is a finite-time screen,
        not a proof of an invariant attracting set or a test of all possible basins.
        """
        parameters = np.asarray(parameters, dtype=self.dtype)
        seeds = np.asarray(seeds, dtype=self.dtype)
        if parameters.ndim != 2 or parameters.shape[1] != 4 or not len(parameters):
            raise ValueError("parameters must have shape (N,4)")
        if seeds.ndim != 2 or seeds.shape[1] != 3 or not len(seeds):
            raise ValueError("seeds must have shape (S,3)")
        if not np.isfinite(parameters).all() or not np.isfinite(seeds).all():
            raise ValueError("parameters and seeds must be finite")
        if not np.isfinite([dt, transient, duration, threshold, escape]).all():
            raise ValueError("integration settings must be finite")
        if dt<=0 or transient<0 or duration<=0 or threshold<=0 or escape<=0:
            raise ValueError("invalid integration settings")
        renorm = 10
        if chunk_steps < renorm or chunk_steps % renorm:
            raise ValueError("chunk_steps must be a positive multiple of 10")
        counts = [int(round(transient/dt)), int(round(duration/(4*dt)))]
        if any(k % renorm for k in counts):
            raise ValueError("transient/dt and duration/(4*dt) must be multiples of 10")
        if not np.isclose(counts[0]*dt, transient) or not np.isclose(counts[1]*dt*4, duration):
            raise ValueError("integration times must be exact multiples of time steps")
        nparam, nseed = len(parameters), len(seeds)
        if initial_states is not None:
            initial_states=np.asarray(initial_states,dtype=self.dtype)
            if (initial_states.ndim!=3 or initial_states.shape[0]!=nparam or
                    initial_states.shape[2]!=3 or initial_states.shape[1]<1 or not np.isfinite(initial_states).all()):
                raise ValueError("initial_states must have finite shape (N,S,3)")
            nseed=initial_states.shape[1]
        n = nparam*nseed
        p = np.ascontiguousarray(np.repeat(parameters,nseed,axis=0))
        state = np.zeros((n,8),dtype=self.dtype)
        state[:,:3] = (np.tile(seeds,(nparam,1)) if initial_states is None
                      else initial_states.reshape(n,3))
        state[:,3:6] = np.asarray([1,2,3])/np.sqrt(14)
        logs = np.zeros((n,4),dtype=self.dtype)
        stats = np.zeros((n,8),dtype=self.dtype)
        stats[:,[0,6]] = np.inf
        stats[:,[1,7]] = -np.inf
        status = np.zeros(n,dtype=np.int32)
        bs,bp,bl,bt,be = [self.buffer(x,readonly=i==1)
                         for i,x in enumerate((state,p,logs,stats,status))]
        total = counts[0]+4*counts[1]
        done=0
        gpu_seconds=0.0
        max_kernel_seconds=0.0
        start=time.perf_counter()
        for window,steps in [(-1,counts[0])]+[(i,counts[1]) for i in range(4)]:
            for offset in range(0,steps,chunk_steps):
                count=min(chunk_steps,steps-offset)
                event=self.advance(self.queue,(n,),None,bs,bp,bl,bt,be,np.int32(n),
                                   np.int32(count),self.scalar(dt),np.int32(window),
                                   np.int32(renorm),self.scalar(escape))
                event.wait()  # Keep Windows GPU submissions short and report real completion.
                elapsed=(event.profile.end-event.profile.start)*1e-9
                gpu_seconds+=elapsed
                max_kernel_seconds=max(max_kernel_seconds,elapsed)
                done+=count
                if progress:
                    progress(done/total,time.perf_counter()-start)
        for array,buf in ((state,bs),(logs,bl),(stats,bt),(status,be)):
            cl.enqueue_copy(self.queue,array,buf).wait()
        wall_seconds=time.perf_counter()-start
        windows=logs.astype(np.float64)/(duration/4)
        lle=windows.mean(axis=1)
        valid=status==0
        windows[~valid]=np.nan
        lle[~valid]=np.nan
        late_range=stats[:,7]-stats[:,6]
        classes=np.full(n,3,dtype=np.uint8)
        classes[(valid)&(lle<=threshold)&(windows[:,-1]<=threshold)]=2
        classes[(valid)&(late_range<1e-5)&(lle<0)]=1
        candidates=valid & (lle>threshold) & (windows.min(axis=1)>threshold/2) & (late_range>1e-4)
        classes[candidates]=4
        classes[status==1]=0
        classes[status==2]=5
        # A double scroll requires repeated visitation of both outer regions.
        # Extrema only provide a weaker descriptive flag, not a topology proof.
        visits_both=(stats[:,0]<-1)&(stats[:,1]>1)&valid
        def shape(x):
            return x.reshape((nparam,nseed)+x.shape[1:])
        return {"lle":shape(lle), "lle_windows":shape(windows),
                "classes":shape(classes), "status":shape(status),
                "visits_both_outer_regions":shape(visits_both),
                "final_states":shape(state[:,:3]), "stats":shape(stats),
                "timing":{"wall_seconds":wall_seconds,"gpu_seconds":gpu_seconds,
                          "max_kernel_seconds":max_kernel_seconds},
                "settings":{"dt":dt,"transient":transient,"duration":duration,
                            "threshold":threshold,"escape":escape,
                            "seeds":seeds.tolist() if initial_states is None else None,
                            "initial_state_mode":"shared" if initial_states is None else "per_parameter",
                            "renormalize_every_steps":renorm,"chunk_steps":chunk_steps}}

    def trace(self, parameters, initial, dt=0.01, steps=1000):
        p=np.atleast_2d(np.asarray(parameters,dtype=self.dtype))
        initial=np.atleast_2d(np.asarray(initial,dtype=self.dtype))
        if p.shape[1]!=4 or initial.shape!=(len(p),3) or steps<1 or dt<=0:
            raise ValueError("invalid trace arguments")
        s=np.zeros((len(p),8),dtype=self.dtype)
        s[:,:3]=initial
        output=np.empty((steps,len(p),4),dtype=self.dtype)
        bs,bp=self.buffer(s),self.buffer(p,True)
        bo=cl.Buffer(self.context,cl.mem_flags.WRITE_ONLY,output.nbytes)
        self.trace_kernel(self.queue,(len(p),),None,bs,bp,bo,np.int32(len(p)),
                          np.int32(steps),self.scalar(dt)).wait()
        cl.enqueue_copy(self.queue,output,bo).wait()
        return output[:,:,:3]
