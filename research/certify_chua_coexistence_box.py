"""Uniform four-cycle persistence by interval multiple shooting.

Each flight has independent starting y,z and duration variables. This avoids
the wrapping of an event-driven single-shooting enclosure over a parameter box.
The Krawczyk operator uses one common center and inverse for all parameters.
"""
import argparse
import json
import time
import numpy as np
from flint import arb, arb_mat, ctx
from certify_chua_coexisting_cycles import ReturnMap, certify as certify_nominal
from chua_interval_map import col, eye, midpoint, bounds
from chua_research_followup import ROOT, OUT, dump

DIR=OUT/'coexistence-box'


def enclosure(record,inward=False):
    """Use endpoint records; the short display ball may be heavily rounded."""
    lo=arb(record['lower']);hi=arb(record['upper'])
    return lo.upper().union(hi.lower()) if inward else lo.lower().union(hi.upper())


def nominal():
    return [arb(453)/50,arb(71729)/5000,-arb(8)/7+arb(3)/500,-arb(5)/7,arb(0)]


class MultipleShooting:
    def __init__(self,record,digits):
        self.original=ReturnMap(record,digits)
        self.schedule=self.original.schedule_rows
        self.n=len(self.schedule)
        self.start_x=[0]+[row[1] for row in self.schedule[:-1]]
        state=col([0,*[arb(str(x)) for x in record['cycle']['q'][:2]]])
        values=[]
        for reg,target,t,sign in self.schedule:
            values.extend([state[1,0],state[2,0],arb(str(t))])
            A,c=self.matrices(nominal(),reg)
            state,_=self.original.flow(A,c,state,arb(str(t)))
            state[0,0]=target
        center=col(values)
        for _ in range(7):
            f,J,_=self.evaluate(center,nominal())
            center=midpoint(center-J.inv()*f)
        self.center=center

    @staticmethod
    def matrices(par,reg):
        alpha,beta,m0,m1,gamma=par
        m=m0 if reg==0 else m1
        A=arb_mat([[-alpha*(1+m),alpha,0],[1,-1,1],[0,-beta,-gamma]])
        b=col([-alpha*(m0-m1)*reg,0,0])
        c=-(A.inv()*b)
        return A,c

    def evaluate(self,q,par):
        n=self.n;F=col([0]*(3*n));J=arb_mat(3*n,3*n);arcs=[]
        for k,(reg,target,guess,sign) in enumerate(self.schedule):
            nxt=(k+1)%n
            start=col([self.start_x[k],q[3*k,0],q[3*k+1,0]])
            T=q[3*k+2,0]
            A,c=self.matrices(par,reg)
            end,E=self.original.flow(A,c,start,T)
            v=A*(end-c)
            dest=col([target,q[3*nxt,0],q[3*nxt+1,0]])
            for i in range(3):
                F[3*k+i,0]=end[i,0]-dest[i,0]
                J[3*k+i,3*k]=E[i,1]
                J[3*k+i,3*k+1]=E[i,2]
                J[3*k+i,3*k+2]=v[i,0]
            J[3*k+1,3*nxt]-=1
            J[3*k+2,3*nxt+1]-=1
            arcs.append(dict(A=A,c=c,start=start,end=end,E=E,v=v,duration=T,region=reg,target=target,sign=sign))
        return F,J,arcs

    def parameter_jacobian(self,q,par):
        """Frechet derivatives of affine exponentials, with Arb enclosures."""
        B=arb_mat(3*self.n,5)
        alpha,beta,m0,m1,gamma=par
        for k,(reg,target,guess,sign) in enumerate(self.schedule):
            m=m0 if reg==0 else m1
            A,c=self.matrices(par,reg)
            H=arb_mat(4,4)
            for i in range(3):
                for j in range(3):H[i,j]=A[i,j]
            H[0,3]=-alpha*(m0-m1)*reg
            derivatives=[arb_mat(4,4) for _ in par]
            derivatives[0][0,0]=-(1+m);derivatives[0][0,1]=1
            derivatives[0][0,3]=-(m0-m1)*reg
            derivatives[1][2,1]=-1
            derivatives[2][0,0]=-alpha if reg==0 else 0
            derivatives[2][0,3]=-alpha*reg
            derivatives[3][0,0]=-alpha if reg!=0 else 0
            derivatives[3][0,3]=alpha*reg
            derivatives[4][2,2]=-1
            initial=col([self.start_x[k],q[3*k,0],q[3*k+1,0],1])
            for p,Hp in enumerate(derivatives):
                block=arb_mat(8,8)
                for i in range(4):
                    for j in range(4):
                        block[i,j]=block[i+4,j+4]=H[i,j]
                        block[i,j+4]=Hp[i,j]
                EE=(block*q[3*k+2,0]).exp()
                dE=arb_mat([[EE[i,j+4] for j in range(4)] for i in range(4)])
                v=dE*initial
                for i in range(3):B[3*k+i,p]=v[i,0]
        return B

    def certify(self,par,radius,verify_path=True,shift=False):
        center=midpoint(self.center)
        center_par=[v.mid() for v in par]
        F0,J0,_=self.evaluate(center,center_par)
        C=midpoint(J0.inv())
        B=self.parameter_jacobian(center,par)
        forcing=-C*F0-(C*B)*col([v-v.mid() for v in par])
        radii=[max(arb('1e-30'),2*v.abs_upper()) for v in forcing.entries()] if radius=='auto' else [arb(radius)]*center.nrows()
        if shift:
            # Keep the original radii: the control must move the proposed box,
            # not enlarge it until it contains the actual nominal orbit again.
            center[0,0]+=arb('0.01')
            F0,J0,_=self.evaluate(center,center_par)
            C=midpoint(J0.inv())
            B=self.parameter_jacobian(center,par)
        box=col([arb(v,r) for v,r in zip(center.entries(),radii)])
        _,J,arcs=self.evaluate(box,par)
        R=eye(J.nrows())-C*J
        K=center-C*F0-(C*B)*col([v-v.mid() for v in par])+R*(box-center)
        contraction=max((sum(abs(R[i,j])*radii[j] for j in range(J.ncols()))/radii[i]).abs_upper() for i in range(J.nrows()))
        inclusion=all(box[i,0].contains_interior(K[i,0]) for i in range(J.nrows()))
        if not(inclusion and contraction<1):
            displacement=max(float(abs(v).abs_upper()) for v in (K-center).entries())
            raise ValueError(f'Uniform Krawczyk inclusion rejected: contraction={float(contraction):.6g}, image radius={displacement:.6g}')
        initial_box=box
        initial_K=K
        for _ in range(4):
            box=col([box[i,0].intersection(K[i,0]) for i in range(J.nrows())])
            _,J,arcs=self.evaluate(box,par)
            R=eye(J.nrows())-C*J
            K=center-C*F0-(C*B)*col([v-v.mid() for v in par])+R*(box-center)
        # Enclose the true return derivative using its independent flight
        # states and times in the validated multiple-shooting box.
        D=eye(2);cells=[];speeds=[];xmax=arb(0);events=[];logdet=arb(0)
        for a in arcs:
            E,v=a['E'],a['v'];assert not v[0,0].contains(0)
            event=arb_mat([[E[i,j]-v[i,0]*E[0,j]/v[0,0] for j in (1,2)] for i in (1,2)])
            D=event*D;events.append(event)
            logdet+=a['A'].trace()*a['duration']
            speeds.append(bounds(v[0,0]));assert a['duration']>0
            if verify_path:
                normalized=dict(start=a['sign']*a['start'],duration=a['duration'],
                    region=abs(a['region']),target=a['sign']*a['target'])
                cells.append(self.original.verify_arc(a['A'],a['sign']*a['c'],normalized)['accepted_cells'])
                for j in range(64):
                    t=arb(j)/64;tt=t.union(arb(j+1)/64)*a['duration']
                    x=self.original.flow(a['A'],a['c'],a['start'],tt)[0][0,0]
                    xmax=max(xmax,x.abs_upper())
        zero=[(i,a) for i,a in enumerate(arcs) if a['target']==0]
        assert len(zero)==2 and zero[0][1]['v'][0,0]<0 and zero[1][1]['v'][0,0]>0
        tr,det=D.trace(),logdet.exp();schur=[1-tr+det,1+tr+det,1-det]
        # Exclude +/-1 without multiplying poorly conditioned event matrices.
        # With det(P)<1, a continuous real 2x2 family can leave the unit disk
        # only through +1 or -1. The nominal cycles were separately certified.
        exclusions=[]
        for multiplier in [-1,1]:
            H=arb_mat(2*self.n,2*self.n)
            for k,event in enumerate(events):
                nxt=(k+1)%self.n
                for i in range(2):
                    for j in range(2):H[2*k+i,2*k+j]=-event[i,j]
                    H[2*k+i,2*nxt+i]+=multiplier if nxt==0 else 1
            pre=midpoint(midpoint(H).inv())
            defect=eye(H.nrows())-pre*H
            norm=max(sum(abs(defect[i,j]) for j in range(H.ncols())).abs_upper() for i in range(H.nrows()))
            if not norm<1:raise ValueError(f'Multiplier {multiplier} exclusion rejected: {float(norm):.6g}')
            exclusions.append(dict(multiplier=multiplier,inverse_defect=bounds(norm)))
        assert det<1
        return dict(passed=True,section_box=[bounds(box[i,0]) for i in [0,1]],
            shooting_box=[bounds(v) for v in initial_box.entries()],krawczyk_image=[bounds(v) for v in initial_K.entries()],
            refined_shooting_box=[bounds(v) for v in box.entries()],
            contraction=bounds(contraction),schur_margins=[bounds(v) for v in schur],
            stability_method='Nominal Schur stability, uniform determinant below one and exclusion of +/-1 by cyclic variational equations.',
            multiplier_exclusions=exclusions,
            trace=bounds(tr),determinant=bounds(det),return_derivative=[[bounds(D[i,j]) for j in range(2)] for i in range(2)],
            period=bounds(sum(a['duration'] for a in arcs)),path_cells=cells,crossing_speeds=speeds,
            x_abs_upper=bounds(xmax),first_upward_return=True)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--digits',type=int,default=80)
    ap.add_argument('--radius',default='1e-8');ap.add_argument('--state-radius',default='1e-4')
    ap.add_argument('--output',default='certificate.json');args=ap.parse_args();ctx.dps=args.digits
    start=time.monotonic();p=nominal();r=arb(args.radius)
    par=[arb(v,r) for v in p[:4]]+[arb(0).union(r)]
    records=json.loads((ROOT/'results/chua-layer-atlas/four-cycles.json').read_text())['records']
    results=[];failure=None
    for row in records:
        f=MultipleShooting(row,args.digits)
        try:
            cert=f.certify(par,args.state_radius)
            base=certify_nominal(row,args.digits)
            _,_,base_arcs=f.original.evaluate([arb(v['ball']) for v in base['root']])
            base_values=[]
            for arc in base_arcs:base_values.extend([arc['start'][1,0],arc['start'][2,0],arc['duration']])
            assert all(enclosure(box,True).contains(value) for box,value in zip(cert['shooting_box'],base_values))
            cert['nominal_stability_certificate']=base
            cert['nominal_orbit_in_shooting_box']=True
            cert['id']=row['id'];results.append(cert)
            print(row['id'],'certified',cert['contraction']['ball'],flush=True)
        except (ValueError,AssertionError,ZeroDivisionError) as error:
            failure=dict(id=row['id'],type=type(error).__name__,message=str(error));break
    gaps=[]
    for i,a in enumerate(results):
        for b in results[i+1:]:
            gap=abs(arb(a['section_box'][0]['ball'])-arb(b['section_box'][0]['ball']))
            assert gap>0
            gaps.append(dict(pair=[a['id'],b['id']],y_gap=bounds(gap)))
    rejected=False
    try:MultipleShooting(records[0],args.digits).certify(par,args.state_radius,False,True)
    except (ValueError,AssertionError,ZeroDivisionError):rejected=True
    result=dict(working_digits=args.digits,parameter_radius=args.radius,state_radius=args.state_radius,
        exact_center=['453/50','71729/5000','-8/7+3/500','-5/7','0'],
        parameter_box={k:bounds(v) for k,v in zip(['alpha','beta','m0','m1','gamma'],par)},
        records=results,separations=gaps,failure=failure,shifted_root_rejected=rejected,
        passed=len(results)==4 and rejected,elapsed_seconds=time.monotonic()-start,
        scope='At least four distinct attracting periodic orbits for every parameter vector in the displayed box; no maximum count, global basin, or hiddenness claim.')
    dump(DIR/args.output,result)
    if not result['passed']:raise SystemExit(str(failure))


if __name__=='__main__':main()
