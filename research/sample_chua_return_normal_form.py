"""Sample a physical return along the secant of two fixed-parameter cycles.

Each event time is solved again for the perturbed initial point. The predicted
scalar curve s -> (s+s^2)/2 has no fitted coefficients after using the two cycle
roots to define s. Results are numerical event-map evaluations.
"""
import argparse
import json
import time
import mpmath as mp
from chua_precision_flights import PrecisionFlights
from study_chua_return_geometry import ROOT, OUT, read


def event_return(row,initial,digits):
    f=PrecisionFlights(row['regions'],row['targets'],True,digits)
    beta=mp.mpf(row['beta']);state=mp.matrix([0,*initial]);arcs=[]
    speeds=[];min_duration=mp.inf;violation=mp.mpf(0)
    for k,(reg,guess) in enumerate(zip(row['regions'],row['q'][2:])):
        slope=f.m0 if reg==0 else f.m1
        A=mp.matrix([[-f.alpha*(1+slope),f.alpha,0],[1,-1,1],[0,-beta,0]])
        c=mp.matrix([reg*(f.m1-f.m0)/(1+f.m1),0,-reg*(f.m1-f.m0)/(1+f.m1)])
        old=state.copy();target=f.targets[k] if k<len(f.targets) else 0
        def at(t):return c+mp.expm(A*t)*(old-c)
        duration=mp.findroot(lambda t:at(t)[0]-target,mp.mpf(guess),solver='newton',
                             df=lambda t:(A*(at(t)-c))[0],tol=mp.mpf(10)**(-(digits-15)),maxsteps=16)
        assert duration>0
        state=at(duration);speed=abs((A*(state-c))[0])
        E=mp.expm(A*duration/32);sample=old.copy()
        for j in range(33):
            violation=max(violation,abs(sample[0])-1 if reg==0 else 1-reg*sample[0])
            sample=c+E*(sample-c)
        speeds.append(speed);min_duration=min(min_duration,duration)
        arcs.append(dict(region=reg,duration=mp.nstr(duration,25),time_change=float(duration-mp.mpf(guess))))
    return -state[1:,:],dict(minimum_crossing_speed=float(min(speeds)),
                             minimum_duration=float(min_duration),
                             maximum_sampled_region_violation=float(violation),
                             maximum_flight_time_change=max(abs(a['time_change']) for a in arcs))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ids',nargs='+',default=['P55','P85'])
    p.add_argument('--digits',type=int,default=75)
    args=p.parse_args()
    for id in args.ids:
        target=OUT/f'{id}-normal-form-d{args.digits}.json'
        if target.exists():print('Retained',target.name,flush=True);continue
        source=ROOT/f'results/chua-knot-atlas/{id}-d70-n192.json'
        if not source.exists():source=OUT/f'{id}-reference.json'
        row=read(source);companion=read(OUT/f'{id}-companion-d70.json')['companion']
        mp.mp.dps=args.digits
        q=mp.matrix([mp.mpf(v) for v in row['q'][:2]])
        d=mp.matrix([mp.mpf(v) for v in companion['q'][:2]])-q
        dd=mp.fdot(d,d);rows=[];started=time.monotonic()
        for s in ['-.5','-.25','0','.25','.5','.75','1','1.25','1.5']:
            value,checks=event_return(row,q+mp.mpf(s)*d,args.digits)
            scalar=mp.fdot(value-q,d)/dd
            off=value-q-scalar*d
            predicted=(mp.mpf(s)+mp.mpf(s)**2)/2
            rows.append(dict(s=float(s),return_coordinate=mp.nstr(scalar,30),
                             predicted=float(predicted),absolute_error=float(abs(scalar-predicted)),
                             transverse_line_error_in_root_separations=float(mp.norm(off)/mp.sqrt(dd)),
                             event_checks=checks))
            print(id,s,'map',mp.nstr(scalar,12),'error',mp.nstr(scalar-predicted,6),flush=True)
        result=dict(id=id,digits=args.digits,rows=rows,elapsed_seconds=time.monotonic()-started,
                    definition='s=0 attracting fixed point; s=1 companion saddle; project the signed half-cycle return along their section secant',
                    predicted_curve='(s+s^2)/2',maximum_absolute_error=max(r['absolute_error'] for r in rows),
                    scope='Local secant samples of the actual event map; not a global conjugacy or an interval normal-form bound')
        target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':main()
