"""Independent ODE checks of selected transverse growth measurements."""
import argparse
import json
import time
import numpy as np
from scipy.integrate import solve_ivp
from trace_chua_window import Chua
from study_chua_return_geometry import OUT, read, norm2, digest


def variational(row, data, step):
    model=Chua()
    beta=float(row['beta'])
    start=row['maximizing_start_index']
    steps=row['maximizing_steps']
    count=len(data['t'])
    times=np.r_[data['t'],data['t']+data['period'],2*data['period']]
    product=np.eye(2)
    errors=[]
    for k in range(start,start+steps):
        i,j=k%count,(k+1)%count
        u=data['xyz'][i]
        B=data['basis'][i]
        def rhs(t,w):
            return np.r_[model.field(beta,w[:3]),
                         (model.jacobian(beta,abs(w[0])<1)@w[3:].reshape(3,2)).ravel()]
        sol=solve_ivp(rhs,[0,times[k+1]-times[k]],np.r_[u,B.ravel()],
                      method='DOP853',rtol=2.5e-14,atol=2.5e-15,max_step=step)
        assert sol.success
        error=np.linalg.norm(sol.y[:3,-1]-data['xyz'][j])
        errors.append(float(error))
        local=data['basis'][j].T@sol.y[3:,-1].reshape(3,2)
        product=local@product
    gain=float(norm2(product))
    return dict(method='DOP853 full field and variation; transverse projection at each mesh point',
                max_step=step,relative_gain_error=abs(gain/row['maximum_gain']-1),
                gain=gain,maximum_local_state_error=max(errors),
                endpoint_matrix=product.tolist(),
                closure_scope='Maximum separate subinterval state error; not a full-period closure error')


def finite_perturbations(row,data,step,amplitudes):
    model=Chua()
    beta=float(row['beta'])
    i=row['maximizing_start_index']
    j=(i+row['maximizing_steps'])%len(data['t'])
    u=data['xyz'][i]
    direction=np.array(row['input_direction'])
    end_B=data['basis'][j]
    outputs=[]
    for eps in amplitudes:
        ends=[]
        for sign in [-1,1]:
            sol=solve_ivp(lambda t,v:model.field(beta,v),[0,row['maximizing_duration']],
                          u+sign*eps*direction,method='DOP853',rtol=2.5e-14,atol=2.5e-15,max_step=step)
            assert sol.success
            ends.append(sol.y[:,-1])
        derivative=end_B.T@(ends[1]-ends[0])/(2*eps)
        gain=float(np.linalg.norm(derivative))
        outputs.append(dict(amplitude=eps,central_difference_gain=gain,
                            relative_error=abs(gain/row['maximum_gain']-1),
                            final_transverse_separation=float(np.linalg.norm(end_B.T@(ends[1]-ends[0])))))
    return dict(method='DOP853 full-field integrations of two perturbed initial states',
                max_step=step,results=outputs,
                scope='Small finite-displacement check over the maximizing interval, not a basin-escape test')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--labels',nargs='+',required=True)
    p.add_argument('--steps',nargs='+',type=float,default=[.025,.0125])
    p.add_argument('--finite',action='store_true')
    args=p.parse_args()
    for label in args.labels:
        path=OUT/(label+'-check.json')
        if path.exists():
            print('Retained',path.name,flush=True)
            continue
        source=OUT/(label+'.json')
        row=read(source)
        with np.load(OUT/(label+'.npz')) as z:
            data={k:z[k] for k in z.files}
        started=time.monotonic()
        report=dict(id=row['id'],source=source.name,source_sha256=digest(source),checks=[])
        for step in args.steps:
            check=variational(row,data,step)
            report['checks'].append(check)
            print(label,step,check,flush=True)
        if args.finite:
            report['finite_perturbations']=finite_perturbations(row,data,args.steps[-1],[1e-8,1e-9,1e-10,1e-11])
            print(label,report['finite_perturbations'],flush=True)
        report['passed']=all(c['relative_gain_error']<2e-5 and c['maximum_local_state_error']<1e-10 for c in report['checks'])
        report['elapsed_seconds']=time.monotonic()-started
        path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
        if not report['passed']:
            raise SystemExit('Independent gain check failed')


if __name__=='__main__':main()
