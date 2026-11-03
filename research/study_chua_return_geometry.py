"""Measure phase-dependent transverse amplification of existing Chua cycles.

The quotient by the flow direction removes infinitesimal phase shifts. Results
use the Euclidean norm of the published dimensionless state, not energy or basin
escape distance. No original periodic root or numerical record is overwritten.
"""
import argparse
import hashlib
import json
import os
import time
from pathlib import Path

import mpmath as mp
import numpy as np

from chua_precision_flights import PrecisionFlights

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get('CHUA_RETURN_OUTPUT', ROOT / 'results/chua-return-geometry')).resolve()
IDS = ['S04', 'S05', 'S06', 'S07', 'P55', 'P65', 'P75', 'P85']


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def digest(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def cross(a, b):
    return mp.matrix([a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]])


def basis(v):
    """Deterministic orthonormal frame in the plane normal to the velocity."""
    tangent = v / mp.sqrt(mp.fdot(v, v))
    axis = mp.zeros(3, 1)
    axis[min(range(3), key=lambda i: abs(tangent[i]))] = 1
    a = cross(tangent, axis)
    a /= mp.sqrt(mp.fdot(a, a))
    b = cross(tangent, a)
    return mp.matrix([[a[i], b[i]] for i in range(3)])


def asarray(m):
    return np.array(m.tolist(), dtype=float)


def mesh(row, subdivisions, digits):
    f = PrecisionFlights(row['regions'], row['targets'], row['symmetric'], digits)
    beta = mp.mpf(row['beta'])
    q = [mp.mpf(v) for v in row['q']]
    ev = f.evaluate(beta, q, False)
    times, points, frames, fields, local = [], [], [], [], []
    elapsed = mp.mpf(0)
    worst_tangent = mp.mpf(0)
    half_product = mp.eye(2)
    for arc in ev['arcs']:
        region = arc['region']
        slope = f.m0 if region == 0 else f.m1
        A = mp.matrix([[-f.alpha*(1+slope), f.alpha, 0], [1, -1, 1], [0, -beta, 0]])
        c = mp.matrix([region*(f.m1-f.m0)/(1+f.m1), 0, -region*(f.m1-f.m0)/(1+f.m1)])
        u = arc['start'].copy()
        dt = arc['duration']/subdivisions
        E = mp.expm(A*dt)
        for k in range(subdivisions):
            v = A*(u-c)
            B = basis(v)
            next_u = c+E*(u-c)
            next_v = A*(next_u-c)
            next_B = basis(next_v)
            L = next_B.T*E*B
            half_product = L*half_product
            worst_tangent = max(worst_tangent, mp.norm(E*v-next_v))
            times.append(float(elapsed+k*dt))
            points.append(list(u))
            frames.append(asarray(B))
            fields.append(list(v))
            local.append(asarray(L))
            u = next_u
        elapsed += arc['duration']
    points = np.array(points, float)
    fields = np.array(fields, float)
    frames = np.array(frames)
    local = np.array(local)
    if row['symmetric']:
        # basis(-v) = basis(v) diag(-1,1) under the convention above.
        D = np.diag([-1., 1.])
        local = np.concatenate([local, D@local@D])
        points = np.concatenate([points, -points])
        fields = np.concatenate([fields, -fields])
        frames = np.concatenate([frames, frames@D])
        times = times + [float(elapsed)+t for t in times]
        dm = mp.matrix([[-1, 0], [0, 1]])
        full_product = dm*half_product*dm*half_product
        period = 2*elapsed
    else:
        full_product = half_product
        period = elapsed
    # The very small second multiplier cannot be obtained by subtracting two
    # nearly equal float64 determinants. Liouville supplies its product instead.
    logdet = ev['logdet']*(2 if row['symmetric'] else 1)
    det = mp.exp(logdet)
    tr = mp.trace(full_product) if hasattr(mp, 'trace') else full_product[0, 0]+full_product[1, 1]
    disc = mp.sqrt(tr*tr-4*det)
    large = (tr+disc)/2 if tr >= 0 else (tr-disc)/2
    expected = mp.mpf('.25') if row['symmetric'] else mp.mpf('-.5')
    # The endpoint mismatch is bounded by the supplied root's closure residual;
    # no frame redefinition is needed at an exact periodic root.
    diag = dict(closure=float(max(abs(v) for v in ev['residual'])),
                dominant_multiplier=mp.nstr(large, 22),
                multiplier_error=float(abs(large-expected)),
                tangent_transport_error=float(worst_tangent),
                determinant=mp.nstr(det, 14), full_product=asarray(full_product).tolist())
    return dict(t=np.array(times), xyz=points, velocity=fields, basis=frames,
                local=local, period=float(period)), diag


def norm2(m):
    a = np.sum(m*m, axis=(-2, -1))
    d = m[..., 0, 0]*m[..., 1, 1]-m[..., 0, 1]*m[..., 1, 0]
    return np.sqrt((a+np.sqrt(np.maximum(a*a-4*d*d, 0)))/2)


def amplification(data):
    local = data['local']
    n = len(local)
    indices = np.arange(n)
    products = np.broadcast_to(np.eye(2), (n, 2, 2)).copy()
    gains = np.ones(n)
    lengths = np.zeros(n, int)
    for k in range(n):
        products = local[(indices+k) % n]@products
        g = norm2(products)
        changed = g > gains
        gains[changed], lengths[changed] = g[changed], k+1
    start = int(np.argmax(gains))
    count = int(lengths[start])
    product = np.eye(2)
    history = [1.]
    maximizing = None
    for k in range(n):
        product = local[(start+k) % n]@product
        history.append(float(norm2(product)))
        if k+1 == count:
            maximizing = product.copy()
    _, singular, vt = np.linalg.svd(maximizing)
    direction = data['basis'][start]@vt[0]
    extended_times = np.r_[data['t'], data['t']+data['period'], 2*data['period']]
    intervals = extended_times[indices+lengths]-extended_times[indices]
    return dict(gain=gains, maximizing_steps=lengths, maximizing_duration=intervals,
                worst_history=np.array(history),
                worst_times=extended_times[start:start+n+1]-extended_times[start]), dict(
        maximum_gain=float(gains[start]), minimum_phase_gain=float(gains.min()),
        phase_gain_ratio=float(gains.max()/gains.min()),
        maximizing_start_index=start, maximizing_steps=count,
        maximizing_phase=float(data['t'][start]/data['period']),
        maximizing_duration=float(intervals[start]),
        input_direction=direction.tolist(), predicted_output_gain=float(singular[0]),
        float64_full_multiplier=float(max(np.linalg.eigvals(products[0]), key=abs).real))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ids', nargs='*')
    parser.add_argument('--subdivisions', type=int, default=8)
    parser.add_argument('--digits', type=int, default=70)
    parser.add_argument('--output', type=Path, default=OUT)
    parser.add_argument('--root-files', nargs='*', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    sources=args.root_files or sorted((ROOT/'results/chua-knot-atlas').glob('*-d70-n192.json'))
    for source in sources:
        source=source.resolve()
        row = read(source)
        if args.ids and row['id'] not in args.ids:
            continue
        label = f"{row['id']}-d{args.digits}-s{args.subdivisions}"
        target = args.output/(label+'.json')
        if target.exists():
            print('Retained', label, flush=True)
            continue
        started = time.monotonic()
        data, diag = mesh(row, args.subdivisions, args.digits)
        curves, summary = amplification(data)
        record = dict(id=row['id'], beta=row['beta'], period=row['period'],
                      symmetric=row['symmetric'], source=source.relative_to(ROOT).as_posix(),
                      source_sha256=digest(source), digits=args.digits,
                      subdivisions_per_affine_flight=args.subdivisions,
                      samples=len(data['t']), diagnostics=diag, **summary,
                      elapsed_seconds=time.monotonic()-started)
        record['passed'] = diag['closure'] < 1e-30 and diag['multiplier_error'] < 1e-20
        np.savez_compressed(args.output/(label+'.npz'), **data, **curves)
        target.write_text(json.dumps(record, indent=2)+'\n', encoding='utf-8', newline='\n')
        print(row['id'], 'G', f"{summary['maximum_gain']:.8g}", 'ratio',
              f"{summary['phase_gain_ratio']:.5g}", 'multiplier error', diag['multiplier_error'],
              'seconds', round(record['elapsed_seconds'], 2), 'passed', record['passed'], flush=True)
        if not record['passed']:
            raise SystemExit('Transverse reconstruction failed its numerical gates.')


if __name__ == '__main__':
    main()
