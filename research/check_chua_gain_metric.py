"""Check gain scaling in the circuit's stored-energy metric."""
import json
import numpy as np
from refine_chua_transient_maxima import frame
from study_chua_return_geometry import OUT, IDS, read, amplification


def main():
    rows=[]
    for id in IDS:
        source=OUT/f'{id}-d70-s8.json';row=read(source)
        with np.load(source.with_suffix('.npz')) as z:
            data={k:z[k] for k in z.files}
        D=np.diag([1/3,1,1/np.sqrt(float(row['beta']))])
        S=np.array([frame(D@v).T@D@B for v,B in zip(data['velocity'],data['basis'])])
        data['local']=np.roll(S,-1,axis=0)@data['local']@np.linalg.inv(S)
        data['basis']=np.array([np.linalg.solve(D,frame(D@v)) for v in data['velocity']])
        curves,summary=amplification(data)
        rows.append(dict(id=id,euclidean_gain=row['maximum_gain'],energy_norm_gain=summary['maximum_gain']))
        print(id,rows[-1],flush=True)
    (OUT/'energy-metric-check.json').write_text(json.dumps(dict(rows=rows,
        norm='sqrt(dx^2/alpha + dy^2 + dz^2/beta)',
        meaning='Amplitude in the quadratic stored-energy metric; the energy ratio is the square of this gain',
        scope='Same coarse phase mesh in both metrics; no continuous maximum or noise model'),indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':main()
