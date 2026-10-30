"""Search off-center transects for families absent from the central census."""
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from chua_gpu import GPU
from detail_analysis import collect_peaks
from survey_chua_layers import recurrence,summarize


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    beta=np.linspace(14.24,14.4,10241);offsets=np.array([-.006,-.003,.003,.006])
    params=np.array([[9.,b,-8/7+d,-5/7] for d in offsets for b in beta]);gpu=GPU();last=[-1]
    def progress(f,s):
        n=int(f*20)
        if n>last[0]:print('scan',n*5,'percent',round(s,1),'seconds',flush=True);last[0]=n
    result=gpu.scan(params,dt=.00125,transient=2000,duration=2400,progress=progress)
    peaks,counts,status=collect_peaks(gpu,params,result['final_states'][:,0],duration=1200,dt=.00125,capacity=512)
    periods=np.array([recurrence(x) if s==0 else 0 for x,s in zip(peaks,status)],dtype=np.int16)
    arrays=dict(parameters=params,beta=beta,m0_offsets=offsets,period=periods,peaks=peaks,peak_counts=counts,peak_status=status,
                **{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
    np.savez_compressed(out/'offcenter-raw.npz',**arrays)
    slices=[]
    for j,d in enumerate(offsets):
        sel=slice(j*len(beta),(j+1)*len(beta))
        slices.append(dict(alpha=9.,m0=float(-8/7+d),**summarize(dict(beta=beta,period=periods[sel],classes=result['classes'][sel],peak_status=status[sel]))))
    report=dict(settings=result['settings'],device=gpu.info,timing=result['timing'],slices=slices)
    (out/'offcenter-summary.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print([(r['m0'],len(r['windows'])) for r in slices],flush=True)


if __name__=='__main__':main()
