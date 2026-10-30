"""Algebra checks of the exact trapped-memory interpretation; no integration."""
import argparse
import json
from pathlib import Path
import numpy as np
from motor_pair import field_jac


def audit():
    rng=np.random.default_rng(2026103087)
    inverse=[];field=[];parent=[];complex_errors=[]
    for _ in range(128):
        rho,r=10**rng.uniform(-2,3,2)
        h,Y,Z=rng.normal(size=3)
        x=(h-r*Y)/rho;c=1-Z;s=-Y
        dh,dY,dZ=field_jac(np.array([h,Y,Z]),np.array([rho,r]))[0]
        v=-rho*x+r*s
        dc=1-c-v*s;ds=-s+v*c
        ref=np.array([(dh-r*dY)/rho,-dZ,-dY])
        field.append(float(np.max(abs(np.array([v,dc,ds])-ref)/(1+abs(ref)))))
        restored=np.array([rho*x-r*s,-s,1-c])
        inverse.append(float(np.max(abs(restored-[h,Y,Z])/(1+abs(restored)))))
        P=r*s;Q=r*(1-c)
        printed=np.array([-P+r*v-v*Q,-Q+v*P])
        transformed=np.array([r*ds,-r*dc])
        parent.append(float(np.max(abs(printed-transformed)/(1+abs(transformed)))))
        zz=c+1j*s
        complex_errors.append(float(abs((dc+1j*ds)-((-1+1j*v)*zz+1))/(1+abs(dc+1j*ds))))
    result=dict(samples=128,seed=2026103087,parameter_sampling='rho,r=10**U(-2,3); standard-normal original motor states',
        inverse_scaled_error=max(inverse),motor_to_memory_field_scaled_error=max(field),
        zero_inertia_parent_field_scaled_error=max(parent),complex_memory_scaled_error=max(complex_errors),
        scope='Random algebraic consistency checks supplement the explicit derivation; no integrations, stability result, or priority claim.')
    assert max(max(inverse),max(field),max(parent),max(complex_errors))<1e-10
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=Path('outputs/motor-memory-algebra.json'))
    args=parser.parse_args();result=audit()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
