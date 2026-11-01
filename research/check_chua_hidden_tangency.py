"""Independent mpmath flights, finite differences and direct-field ODE controls."""
import json
import argparse
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
from flint import arb
from chua_research_followup import OUT,dump
from chua_lobe_map import LobeMap
from certify_chua_hidden_tangency import main as certificate


def main(output_directory=None):
    destination=output_directory or OUT/'hidden-search'
    source=json.loads((OUT/'hidden-search/connection-center-14-2.json').read_text())
    cert=json.loads((OUT/'hidden-search/hidden-tangency-160.json').read_text())
    mp.mp.dps=90;p=source['parameters'];model=LobeMap(p)
    alpha,m0,m1=[mp.mpf(str(p[k])) for k in ['alpha','m0','m1']]
    b0=mp.mpf(source['beta']);a0=mp.mpf(source['a']);V=mp.matrix(source['coordinates']);W=V**-1
    def matrices(beta,reg):
        m=m1 if reg else m0;x=(m1-m0)/(1+m1) if reg else 0
        return mp.matrix([[-alpha*(1+m),alpha,0],[1,-1,1],[0,-beta,0]]),mp.matrix([x,0,-x])
    def apply(q,beta):
        guesses=model.apply(np.array([float(v) for v in q]),float(beta))[2]['arcs']
        state=mp.matrix([1,*q]);M=mp.eye(3);period=mp.mpf(0)
        for reg,guess in zip([1,0],guesses):
            A,c=matrices(beta,reg);t=mp.mpf(str(guess['duration']))
            for _ in range(12):
                E=mp.expm(A*t);end=c+E*(state-c);v=A*(end-c);change=(end[0]-1)/v[0];t-=change
                if abs(change)<mp.mpf('1e-80'):break
            else:raise ValueError('Independent crossing Newton exhausted')
            E=mp.expm(A*t);state=c+E*(state-c);state[0]=1;M=E*M;period+=t
        v=A*(state-c);J=((mp.eye(3)-v*mp.matrix([[1,0,0]])/v[0])*M)[1:3,1:3]
        return state[1:3,0],J,period
    def fixed(beta):
        q=mp.matrix(source['parent_guess'])
        for _ in range(12):
            y,J,_=apply(q,beta);r=y-q
            if mp.norm(r,mp.inf)<mp.mpf('1e-75'):return q
            q-=(J-mp.eye(2))**-1*r
        raise ValueError('Independent saddle Newton exhausted')
    traces=[]
    def gap(a,beta,trace=False):
        saddle=fixed(beta);q=saddle+V*mp.matrix([a,0])
        for k in range(source['iterations']):
            old=q.copy();q,J,T=apply(q,beta)
            if trace:traces.append(dict(index=k,start=[float(v) for v in old],end=[float(v) for v in q],duration=float(T)))
        return W*(q-saddle)
    center=gap(a0,b0,True);print('independent connection evaluated',flush=True)
    tests=[]
    for h in [mp.mpf('1e-15'),mp.mpf('5e-16')]:
        ga=gap(a0+h,b0);gb=gap(a0-h,b0)
        curvature=(ga[0]-2*center[0]+gb[0])/(h*h)
        slope=(ga[0]-gb[0])/(2*h)
        hb=mp.mpf('1e-12');unfold=(gap(a0,b0+hb)[0]-gap(a0,b0-hb)[0])/(2*hb)
        assert arb(cert['normal_curvature']['ball']).contains(arb(str(curvature)))
        assert arb(cert['parameter_unfolding']['ball']).contains(arb(str(unfold)))
        tests.append(dict(step=str(h),curvature=str(curvature),slope=str(slope),unfolding=str(unfold)))
        print('independent differences',tests[-1],flush=True)
    assert abs(center[0])<mp.mpf('1e-25')
    ode=[]
    def field(t,v):
        x,y,z=v;h=float(m1)*x+float(m0-m1)*(abs(x+1)-abs(x-1))/2
        return [float(alpha)*(y-x-h),x-y+z,-float(b0)*y]
    for k in [0,18,35]:
        row=traces[k];sol=solve_ivp(field,[0,row['duration']],[1,*row['start']],method='DOP853',rtol=2e-12,atol=2e-14,max_step=.002)
        error=float(max(abs(sol.y[:,-1]-np.r_[1,row['end']])))
        assert sol.success and error<2e-8
        ode.append(dict(return_index=k,error=error,steps=len(sol.t)))
    # Shift the whole tested beta interval, retaining the original a center.
    negative=[]
    for shift in ['-1e-9','1e-9']:
        target=destination/('shifted-tangency-'+('minus' if shift[0]=='-' else 'plus')+'.json')
        result=certificate(140,shift,'1e-15','5e-12','8e-13',target)
        assert result['status']=='failed'
        dump(target,result);negative.append(dict(shift=shift,failed_gates=[k for k,v in result['gates'].items() if not v]))
    dump(destination/'independent-tangency-checks.json',dict(passed=True,working_digits=90,
        center=[str(v) for v in center],finite_differences=tests,ode=ode,negative_controls=negative,section_trace=traces,
        scope='Independent numerical controls of the enclosure implementation; not a second proof or an external review.'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-directory',type=Path)
    main(parser.parse_args().output_directory)
