"""Sample the long-parent return away from its original attracting root."""
import mpmath as mp
import numpy as np
from test_chua_nested_returns import OUT, read, save
from chua_block_return import BlockShooting
from chua_block_return import BlockReturn
from sample_chua_return_normal_form import event_return


def main():
    rows=[]
    for cal in read(OUT/'calibration.json')['parents']:
        shoot=BlockShooting(cal,1,65);curves=[];controls=[]
        for M in ['1.75','1.95','2']:
            samples=[]
            for X in np.linspace(-2.05,2.05,83):
                x=mp.mpf(str(X));s=-2*x-mp.mpf('.5')
                ev=shoot.block([s,0],mp.mpf(M));actual=-mp.mpf('.25')-ev['normalized'][0]/2
                samples.append(dict(x=float(x),actual=float(actual),quadratic=float(mp.mpf(M)-x*x),
                    off_secant=float(ev['normalized'][1]),time=float(ev['time']),minimum_speed=float(ev['minimum_speed'])))
            curves.append(dict(M=float(M),beta=mp.nstr(shoot.beta(mp.mpf(M)),50),samples=samples,
                maximum_scalar_error=max(abs(s['actual']-s['quadratic']) for s in samples),
                maximum_off_secant=max(abs(s['off_secant']) for s in samples)))
        for M,X in [('1.75','-2'),('1.75','0'),('1.75','2'),('2','-2'),('2','0'),('2','2')]:
            s=-2*mp.mpf(X)-mp.mpf('.5');state=shoot.origin+shoot.S*mp.matrix([s,0])
            row=dict(cal['representative'],beta=mp.nstr(shoot.beta(mp.mpf(M)),60))
            independent,checks=event_return(row,state,80)
            mp.mp.dps=65;ev=shoot.block([s,0],mp.mpf(M))
            controls.append(dict(M=float(M),X=float(X),independent_error=float(mp.norm(ev['point']-independent)),checks=checks))
        children=[]
        for k in [3,5,7]:
            child=read(OUT/f"{cal['id']}-k{k}-reference.json")
            x=np.array([-.25-.5*float(p[0]) for p in child['block_coordinates']]);y=np.array([float(p[1]) for p in child['block_coordinates']]);M=float(child['normalized_coordinate'])
            children.append(dict(id=child['id'],M=M,x=x.tolist(),off_secant=y.tolist(),
                maximum_quadratic_residual=float(max(abs(np.roll(x,-1)-(M-x*x))))))
        result=dict(parent=cal['id'],curves=curves,independent_controls=controls,children=children)
        rows.append(result);print(cal['id'],'map errors',[c['maximum_scalar_error'] for c in curves],
            'off secant',[c['maximum_off_secant'] for c in curves],flush=True)
    save(OUT/'parent-block-map.json',dict(rows=rows,
        coordinates='X=-1/4-s/2, with s along the original attracting-to-companion secant; M=(beta_fold-beta)/parent_width-1/4',
        scope='Finite secant samples of a two-dimensional return. No fitted correction, invariant-curve theorem, or global conjugacy.'))


def moving_coordinates():
    """Post-search diagnostic; it must never replace the frozen forecasts."""
    rows=[]
    for parent in ['S06','S11-n3']:
        g=read(OUT/f'{parent}-k7-geometry.json');r=g['parent'];other=g['continued_attracting_branch']
        f=BlockReturn(r,70);beta=mp.mpf(r['beta']);q=mp.matrix([mp.mpf(v) for v in r['q'][:2]])
        d=mp.matrix([mp.mpf(v) for v in other['q'][:2]])-q;dd=mp.fdot(d,d)
        ev=f.evaluate(beta,q);tr=ev['matrix'][0,0]+ev['matrix'][1,1]
        mu=(tr+mp.sqrt(tr*tr-4*ev['determinant']))/2
        lo=-mu/2;hi=mu/2-1;M=lo+lo*lo;samples=[]
        for X in np.linspace(-2.1,2.1,85):
            x=mp.mpf(str(X));s=(x-lo)/(hi-lo);e=f.evaluate(beta,q+s*d)
            scalar=mp.fdot(e['point']-q,d)/dd;actual=lo+scalar*(hi-lo)
            off=mp.norm(e['point']-q-scalar*d)/mp.sqrt(dd)
            samples.append(dict(x=float(x),actual=float(actual),quadratic=float(M-x*x),off_secant=float(off)))
        rows.append(dict(parent=parent,beta=r['beta'],parent_half_multiplier=mp.nstr(mu,35),M=float(M),samples=samples,
            maximum_scalar_error=max(abs(s['actual']-s['quadratic']) for s in samples),
            maximum_off_secant=max(s['off_secant'] for s in samples)))
        print(parent,'moving-frame error',rows[-1]['maximum_scalar_error'],flush=True)
    save(OUT/'moving-coordinate-map.json',dict(rows=rows,
        scope='Post-search diagnostic using the two parent fixed points and the companion multiplier at each k=7 beta. Not a prospective forecast or a conjugacy proof.',
        coordinates='X_C=-mu_C/2, X_A=mu_C/2-1, M=X_C+X_C^2; linear interpolation along the same-beta C-to-A secant'))


if __name__=='__main__':
    import sys
    if '--moving' in sys.argv:moving_coordinates()
    else:main()
