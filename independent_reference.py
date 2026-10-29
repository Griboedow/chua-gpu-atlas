"""Independent DOP853 reference using the full piecewise vector field.

Every Runge-Kutta stage chooses its actual diode branch. Nonterminal crossing
events are diagnostic: missing a grazing pair cannot leave the dynamics on
an extrapolated affine branch. Tangent normalization uses a scaled norm.
"""
import time
import numpy as np
from scipy.integrate import solve_ivp
from chua_gpu import rhs

def reference_lle(parameters,initial,transient=40,duration=160):
    if transient<0 or duration<=0 or not np.isfinite([transient,duration]).all():
        raise ValueError('Invalid observation times')
    s=np.r_[np.asarray(initial,dtype=float),np.array([1,2,3])/np.sqrt(14)]
    logs=np.zeros(4);maxabs=0.;crossings=0;started=time.perf_counter();current=0.
    def minus(t,v):return v[0]+1
    def plus(t,v):return v[0]-1
    # Events do not change or select the vector field.
    events=[minus,plus]
    for window,length in [(-1,transient)]+[(i,duration/4) for i in range(4)]:
        end=current+length
        while current<end-1e-12:
            target=min(end,current+.25)
            sol=solve_ivp(lambda t,v:rhs(t,v,parameters),(current,target),s,
                method='DOP853',rtol=2e-9,atol=2e-11,
                max_step=min(.1,.5/np.sqrt(parameters[1])),events=events)
            if not sol.success:raise RuntimeError(sol.message)
            s=sol.y[:,-1];current=target
            maxabs=max(maxabs,float(np.max(abs(sol.y[:3]))))
            if not np.isfinite(s).all() or maxabs>1e6:
                raise RuntimeError('Independent trajectory escaped the 1e6 observation cutoff')
            crossings+=sum(len(times) for times in sol.t_events)
            scale=float(np.max(abs(s[3:])))
            if not scale>0:raise RuntimeError('Degenerate tangent vector')
            normalized=s[3:]/scale;unitnorm=np.linalg.norm(normalized)
            if window>=0:logs[window]+=np.log(scale)+np.log(unitnorm)
            s[3:]=normalized/unitnorm
        print(f'Full-field DOP853: t={current:g}',flush=True)
    return {'lle':float(logs.sum()/duration),'windows':(logs/(duration/4)).tolist(),
        'recorded_crossings':crossings,'max_absolute_coordinate':maxabs,
        'transient':transient,'duration':duration,'wall_seconds':time.perf_counter()-started,
        'rtol':2e-9,'atol':2e-11,'renormalization_interval':.25,
        'method':'DOP853 evaluates the full piecewise state and tangent field at every stage; nonterminal x=±1 event diagnostics.',
        'event_scope':'Recorded events are diagnostic, not a guarantee of detecting every grazing crossing pair.'}
