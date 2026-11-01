"""Independent full-field checks of the newly enclosed periodic orbits."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from certify_chua_coexisting_cycles import OUT


def numeric(x):
    return float(x['midpoint'])


def check_cycle(cycle, parameters):
    alpha,beta,m0,m1=parameters
    u=np.r_[0.,[numeric(v) for v in cycle['root']]]
    period=numeric(cycle['period'])
    def field(t,w):
        x,y,z=w[:3]
        h=m1*x+(m0-m1)*(abs(x+1)-abs(x-1))/2
        f=np.array([alpha*(y-x-h),x-y+z,-beta*y])
        if len(w)==3:return f
        slope=m0 if abs(x)<1 else m1
        A=np.array([[-alpha*(1+slope),alpha,0],[1,-1,1],[0,-beta,0]])
        return np.r_[f,(A@w[3:].reshape(3,3)).ravel()]
    checks=[]
    for method,step in [('DOP853',.006),('Radau',.012)]:
        sol=solve_ivp(field,[0,period],u,method=method,rtol=2e-12,atol=2e-13,max_step=step)
        closure=float(np.max(abs(sol.y[:,-1]-u)))
        assert sol.success and closure<2e-8,(cycle['id'],method,closure)
        checks.append(dict(method=method,step=step,closure=closure))
    var=solve_ivp(field,[0,period],np.r_[u,np.eye(3).ravel()],method='DOP853',rtol=2e-11,atol=2e-12,max_step=.004)
    M=var.y[3:,-1].reshape(3,3); v=field(0,var.y[:3,-1])
    J=((np.eye(3)-np.outer(v,[1,0,0])/v[0])@M)[1:,1:]
    rho=float(max(abs(np.linalg.eigvals(J))))
    enclosed=max(abs(numeric(x)) for x in cycle['multipliers'])
    error=abs(rho-enclosed)
    assert var.success and rho<1 and error<2e-6
    return dict(id=cycle['id'],checks=checks,variational_radius=rho,interval_radius_midpoint=enclosed,
                multiplier_error=error,passed=True)


def main():
    four=json.loads((OUT/'four-cycle-certificate-80.json').read_text())
    hidden=json.loads((OUT/'hidden-cycle-certificate-80.json').read_text())
    stability=json.loads((OUT/'stability-certificate-80.json').read_text())
    records=[]
    for c in four['records']:
        records.append(check_cycle(c,[9.06,14.3458,-8/7+.006,-5/7]))
        print(c['id'],'independent checks passed',flush=True)
    for c,pars in [(hidden['periodic_orbit'],[8.4,12,-1.2,-.05]),
                   (stability['unique_equilibrium_hidden_cycle']['periodic_orbit'],[8.4,12,-.95,-.05])]:
        records.append(check_cycle(c,pars)); print(c['id'],'independent checks passed',flush=True)
    # Compare repeated certificates using numeric INTERVALS, not text equality.
    from flint import arb
    repeat=[]
    for name in ['four-cycle-certificate','hidden-cycle-certificate','stability-certificate']:
        a=json.loads((OUT/(name+'-80.json')).read_text())
        b=json.loads((OUT/(name+'-100.json')).read_text())
        def recurse(x,y,path=''):
            if isinstance(x,dict) and 'midpoint' in x and 'radius' in x:
                # `ball` includes serialization rounding. Midpoint/radius are
                # diagnostic strings, not a lossless serialized interval pair.
                # Computed upper/lower bounds need not agree across precisions.
                if path.endswith('/contraction'):
                    assert arb(x['ball'])<1 and arb(y['ball'])<1
                elif path.endswith('/periodic_x_amplitude_upper'):
                    assert arb(x['ball'])<12 and arb(y['ball'])<12
                elif 'excluded_ball_radii' in path or path.endswith('/excluded_origin_radius'):
                    assert arb(x['ball'])>0 and arb(y['ball'])>0
                else:
                    assert arb(x['ball']).overlaps(arb(y['ball'])),path
            elif isinstance(x,dict):
                for k in x:
                    if k!='working_digits':recurse(x[k],y[k],path+'/'+k)
            elif isinstance(x,list):
                assert len(x)==len(y)
                for i,(xx,yy) in enumerate(zip(x,y)):recurse(xx,yy,path+'/'+str(i))
        recurse(a,b); repeat.append(dict(certificate=name,intervals_overlap=True))
    result=dict(passed=True,orbits=records,precision_repeats=repeat,
                scope='Independent numerical controls support the interval implementation; they do not replace its analytical justification.')
    (OUT/'independent-checks.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()
