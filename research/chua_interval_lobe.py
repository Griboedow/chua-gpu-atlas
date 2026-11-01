"""Arb return jets for two transverse flights on the x=1 lobe section."""
import numpy as np
from flint import arb,arb_mat,ctx
from chua_interval_map import IntervalHalfMap,col,eye,midpoint
from chua_lobe_map import LobeMap
from certify_chua_layer_organizer import submatrix


class IntervalLobe(IntervalHalfMap):
    moving_frame=True
    regions=(1,0)
    targets=(1,1)
    section_x=1
    return_sign=1

    def __init__(self,parameters,digits=80):
        ctx.dps=digits;self.parameters=parameters
        self.alpha,self.m0,self.m1=[arb(str(parameters[k])) for k in ['alpha','m0','m1']]
        self.float_map=LobeMap(parameters)

    def matrices(self,beta,region):
        slope=self.m0 if region==0 else self.m1
        A=arb_mat([[-self.alpha*(1+slope),self.alpha,0],[1,-1,1],[0,-beta,0]])
        x=region*(self.m1-self.m0)/(1+self.m1)
        return A,col([x,0,-x]),arb_mat([[0,0,0],[0,0,0],[0,-1,0]])

    def schedule(self,point,beta):
        state=np.array([1.,float(point[0].mid()),float(point[1].mid())]);times=[]
        mats=self.float_map.matrices(float(beta.mid()))
        for reg in self.regions:
            dt,target=self.float_map.event(state,reg,mats[reg])
            if target!=1:raise ValueError('Lobe itinerary left x>-1')
            A,c,e,V,W=mats[reg];state=c+((V*np.exp(e*dt))@W).real@(state-c);state[0]=1
            times.append(dt)
        return times


def saddle_center(f,beta,guess):
    p=col([arb(str(v)) for v in guess])
    for _ in range(12):
        value,J,_,_=f.evaluate(p.entries(),beta)
        if all(abs(v)<arb(10)**(-(ctx.dps-25)) for v in (value-p).entries()):break
        p=midpoint(p-(submatrix(J)-eye(2)).inv()*(value-p))
    return p
