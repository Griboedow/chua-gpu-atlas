"""Independent event-resolved Lyapunov checks for selected far-field candidates."""
import time
import numpy as np
from scipy.integrate import solve_ivp
from explore_extent import OUT
from run import save_json

def event_lle(parameters,initial,transient=40,duration=160):
    a,b,m0,m1=parameters
    s=np.r_[initial,np.array([1,2,3])/np.sqrt(14)]
    region=1 if s[0]>1 else (-1 if s[0]<-1 else 0)
    current=0.;logs=np.zeros(4);crossings=0;maxabs=0.;started=time.perf_counter()
    for block in range(int(transient+duration)):
        end=float(block+1)
        while current<end-1e-11:
            slope=m0 if region==0 else m1
            offset=region*(m0-m1)
            jac=np.array([[-a*(1+slope),a,0],[1,-1,1],[0,-b,0.]])
            forcing=np.array([-a*offset,0,0])
            def flow(t,v):
                return np.r_[jac@v[:3]+forcing,jac@v[3:]]
            if region==0:
                transitions=[(-1.,-1,-1),(1.,1,1)]
            elif region==1:
                transitions=[(1.,-1,0)]
            else:
                transitions=[(-1.,1,0)]
            events=[]
            for boundary,direction,new_region in transitions:
                def event(t,v,boundary=boundary):return v[0]-boundary
                event.direction=direction;event.terminal=True;events.append(event)
            sol=solve_ivp(flow,(current,end),s,method='DOP853',rtol=2e-9,atol=2e-11,
                          max_step=min(.1,.5/np.sqrt(b)),events=events)
            if not sol.success:raise RuntimeError(sol.message)
            previous_time=current;current=float(sol.t[-1]);s=sol.y[:,-1]
            maxabs=max(maxabs,float(np.max(np.abs(sol.y[:3]))))
            event_indices=[i for i,times in enumerate(sol.t_events) if len(times)]
            if event_indices:
                which=event_indices[0];region=transitions[which][2];crossings+=1
            if current<=previous_time+1e-14:
                raise RuntimeError('Zero-time boundary retrigger')
        norm=np.linalg.norm(s[3:])
        if block>=transient:
            logs[min(3,int((block-transient)//(duration/4)))]+=np.log(norm)
        s[3:]/=norm
        if (block+1)%50==0:print(f'Event-resolved reference: t={block+1}',flush=True)
    return {'lle':float(logs.sum()/duration),'windows':(logs/(duration/4)).tolist(),
            'boundary_crossings':crossings,'max_absolute_coordinate':maxabs,
            'transient':transient,'duration':duration,'wall_seconds':time.perf_counter()-started}

if __name__=='__main__':
    data=np.load(OUT/'verified.npz')
    keep=(data['classes']==4).any(axis=1)
    indices=np.flatnonzero(keep)
    chosen=np.unique([indices[np.argmax(data['parameters'][indices,0])],
                      indices[np.argmax(data['parameters'][indices,1])],
                      indices[np.argmin(data['parameters'][indices,2])]])
    results=[]
    for idx in chosen:
        p=data['parameters'][idx]
        seed=int(np.argmax(np.where(data['classes'][idx]==4,data['lle'][idx],-np.inf)))
        print('Checking parameters',p,'seed',seed,flush=True)
        result=event_lle(p,data['final_states'][idx,seed])
        result.update({'parameters':p.tolist(),'gpu_lle':float(data['lle'][idx,seed]),'seed_index':seed})
        results.append(result)
        save_json(OUT/'event_resolved_checks.json',{'method':'DOP853 with exact detection of x=±1 crossings; constant affine field within each region. Tangent vector stays continuous because the vector field is continuous.',
                   'checks':results,'scope':'Selected examples only; starts from the end of the GPU observation, then discards an additional transient.'})
        print(result,flush=True)
