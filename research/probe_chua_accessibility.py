"""Finite-time basin screen near every equilibrium of the four-cycle example.

This is deliberately not a proof of hiddenness or of exhaustive basin coverage.
Classifications use the last upward section crossing, with unresolved outcomes
retained. Step refinement reuses exactly the same initial conditions.
"""
import argparse
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/chua-accessibility"


def field(u, alpha, beta, m0, m1):
    x, y, z = u.T
    return np.column_stack((alpha*(y-(1+m1)*x-(m0-m1)*np.clip(x,-1,1)), x-y+z, -beta*y))


def simulate(initial, parameters, step, duration, roots):
    u = initial.copy()
    section = np.full((len(u), 2), np.nan)
    for k in range(round(duration/step)):
        a = field(u, **parameters)
        b = field(u+step*a/2, **parameters)
        c = field(u+step*b/2, **parameters)
        d = field(u+step*c, **parameters)
        v = u+step*(a+2*b+2*c+d)/6
        crossing = (u[:,0] < 0) & (v[:,0] >= 0)
        w = -u[crossing,0]/(v[crossing,0]-u[crossing,0])
        if np.any(crossing):
            left,right=u[crossing],v[crossing]
            fl=a[crossing]; fr=field(right,**parameters)
            # Cubic Hermite dense interpolation, followed by a scalar Newton
            # solve for x=0. Linear interpolation can misclassify large cycles.
            c0=left; c1=step*fl
            c2=3*(right-left)-step*(2*fl+fr)
            c3=2*(left-right)+step*(fl+fr)
            for _ in range(4):
                value=c0[:,0]+w*(c1[:,0]+w*(c2[:,0]+w*c3[:,0]))
                derivative=c1[:,0]+w*(2*c2[:,0]+3*w*c3[:,0])
                w=np.clip(w-value/derivative,0,1)
            hit=c0+w[:,None]*(c1+w[:,None]*(c2+w[:,None]*c3))
            section[crossing]=hit[:,1:]
        u = v
    distance = np.linalg.norm(section[:,None,:]-roots[None,:,:], axis=2)
    valid = np.isfinite(distance).all(axis=1)
    labels = np.full(len(u), -1)
    labels[valid] = np.argmin(distance[valid], axis=1)
    best = np.full(len(u), np.inf)
    best[valid] = np.min(distance[valid], axis=1)
    # Less than one quarter of the closest separation between the four roots.
    labels[best > 1e-4] = -1
    return labels, best, section, u


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--duration", type=float, default=600)
    p.add_argument("--directions", type=int, default=64)
    args = p.parse_args()
    prior = json.loads((ROOT/"results/chua-layer-atlas/four-cycles.json").read_text())
    parameters = {k: prior[k] for k in ("alpha", "beta", "m0", "m1")}
    roots = np.array([r["cycle"]["q"][:2] for r in prior["records"]])
    outer = (prior["m1"]-prior["m0"])/(1+prior["m1"])
    equilibria = np.array([[0,0,0], [outer,0,-outer], [-outer,0,outer]])
    i = np.arange(args.directions)
    z = 1-2*(i+.5)/args.directions
    phi = np.pi*(3-np.sqrt(5))*i
    directions = np.column_stack((np.sqrt(1-z*z)*np.cos(phi), np.sqrt(1-z*z)*np.sin(phi), z))
    radii = [1e-2, 1e-4, 1e-6, 1e-8]
    initial = np.array([e+r*d for e in equilibria for r in radii for d in directions])
    records, arrays = [], dict(initial=initial, equilibria=equilibria, radii=np.array(radii))
    for step in (.01, .005):
        labels, errors, sections, final = simulate(initial, parameters, step, args.duration, roots)
        counts = []
        cube = labels.reshape(3,len(radii),args.directions)
        for j in range(3):
            for k,r in enumerate(radii):
                counts.append(dict(equilibrium=j,radius=r,counts=[int(np.sum(cube[j,k]==label)) for label in (-1,0,1,2,3)]))
        records.append(dict(step=step, counts=counts))
        tag = str(step)
        arrays.update({tag+"_labels":labels,tag+"_section":sections,tag+"_errors":errors})
        print("step",step,"counts",np.unique(labels,return_counts=True),flush=True)
    different = arrays["0.01_labels"] != arrays["0.005_labels"]
    report = dict(parameters=parameters, duration=args.duration, directions=args.directions,section_interpolation='cubic Hermite',
                  radii=radii, records=records, refinement_disagreements=int(different.sum()),
                  initial_count=len(initial), label_order=["unresolved","F1a","F1b","F2a","F2b"],
                  scope="Finite-radius finite-time numerical screen. Neither hiddenness nor self-excitation is certified.")
    OUT.mkdir(exist_ok=True)
    (OUT/"equilibrium-screen.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    np.savez_compressed(OUT/"equilibrium-screen.npz", **arrays)


if __name__ == "__main__":
    main()
