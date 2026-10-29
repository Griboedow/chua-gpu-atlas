"""Direct-quadrature DOP853 controls, independent of GPU lookup/integrator."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
from quantum_neuron import rhs

def integrate(p,start,transient,duration,trace=False):
    state=np.r_[start,np.array([1,2,3,4])/np.sqrt(30)]
    windows=np.zeros(4);time=0.
    settings=dict(method='DOP853',rtol=2e-10,atol=2e-12,max_step=min(.1,.5/p[3]))
    for k in range(int(transient+duration)):
        sol=solve_ivp(lambda t,y:rhs(t,y,p,1024),(time,time+1),state,**settings)
        if not sol.success:raise RuntimeError(sol.message)
        state=sol.y[:,-1];norm=np.linalg.norm(state[4:]);state[4:]/=norm
        if k>=transient:windows[min(3,int(4*(k-transient)/duration))]+=np.log(norm)
        time+=1
    result={'parameters':p.tolist(),'start':np.asarray(start).tolist(),'transient':transient,'duration':duration,
            'lle_windows':(windows/(duration/4)).tolist(),'final_state':state[:4].tolist(),
            'method':'DOP853; direct 1024-node Gauss-Hermite; normalize every 1 time unit','rtol':settings['rtol'],'atol':settings['atol']}
    if trace:
        t=np.linspace(0,40,8001)
        sol=solve_ivp(lambda t,y:rhs(t,y,p,1024),(0,40),state[:4],t_eval=t,**settings)
        if not sol.success:raise RuntimeError(sol.message)
        peaks,_=find_peaks(sol.y[3]);result['voltage_span']=float(np.ptp(sol.y[3]));result['voltage_peak_span']=float(np.ptp(sol.y[3,peaks]));result['max_bloch_radius']=float(np.sqrt((sol.y[:3]**2).sum(axis=0)).max())
        return result,t,sol.y
    return result

if __name__=='__main__':
    path=Path('results/candidate-research');results=[]
    p=np.array([1.73,5,.1,7,1,1,.8,0.])
    result,t,y=integrate(p,[.05,0,.95,.1],300,200,trace=True)
    print('Paper periodic control',json.dumps(result),flush=True);results.append({'case':'paper_periodic_amplitude',**result})
    np.savez_compressed(path/'periodic_reference.npz',t=t,state=y,parameters=p)
    a=np.load(path/'refined/scan.npz')
    # The two irregular-looking scout histories, now starting from their long-run endpoints.
    for index in [6,9]:
        result,t,y=integrate(a['parameters'][index],a['final_states'][index],100,400,trace=True)
        print('Refined endpoint',index,json.dumps(result),flush=True);results.append({'case':f'refined_endpoint_{index}',**result})
        np.savez_compressed(path/f'reference_{index}.npz',t=t,state=y,parameters=a['parameters'][index])
    (path/'cpu_reference.json').write_text(json.dumps(results,indent=2))
