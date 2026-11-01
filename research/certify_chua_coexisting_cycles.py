"""Validate four first-return sinks at one exact Chua parameter point.

The old numerical records supply guesses only. Arb event inclusions, full-time
itinerary checks, Krawczyk inclusion and the real 2x2 Schur criterion establish
the new claims. This program does not classify their basins as hidden.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq
from flint import arb, arb_mat, ctx
from chua_interval_map import IntervalHalfMap, col, eye, midpoint, bounds

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/chua-accessibility"


class ReturnMap(IntervalHalfMap):
    def __init__(self, record, digits, parameters=None):
        ctx.dps = digits
        if parameters is None:
            self.alpha = arb(453) / 50
            self.beta = arb(71729) / 5000
            self.m0 = -arb(8) / 7 + arb(3) / 500
            self.m1 = -arb(5) / 7
        else:
            self.alpha, self.beta, self.m0, self.m1 = parameters
        self.record = record
        # Split the central flight across zero. Each resulting central flight
        # remains in one signed half-space, which certifies the first return.
        q = record["cycle"]["q"]
        state = np.array([0., *q[:2]])
        self.schedule_rows = []
        for reg, target, duration in zip(record["regions"], [*record["targets"], 0], q[2:]):
            A, c, _ = self.matrices(None, reg)
            An = np.array([[float(A[i, j].mid()) for j in range(3)] for i in range(3)])
            cn = np.array([float(x.mid()) for x in c.entries()])
            start_x = int(round(state[0]))
            if reg == 0 and start_x * target < 0:
                crossing = brentq(lambda t: (cn + expm(An*t)@(state-cn))[0], 0., duration)
                self.schedule_rows.append((reg, 0, crossing, start_x))
                self.schedule_rows.append((reg, target, duration-crossing, target))
            else:
                sign = reg if reg else (target or start_x)
                self.schedule_rows.append((reg, target, duration, sign))
            state = cn + expm(An*duration)@(state-cn)
            state[0] = target

    def matrices(self, unused, reg):
        m = self.m0 if reg == 0 else self.m1
        A = arb_mat([[-self.alpha*(1+m), self.alpha, 0], [1, -1, 1], [0, -self.beta, 0]])
        x = reg*(self.m1-self.m0)/(1+self.m1)
        return A, col([x, 0, -x]), None

    def evaluate(self, point, verify_path=False):
        state = col([0, *point])
        D = arb_mat([[0, 0], [1, 0], [0, 1]])
        arcs = []
        for reg, target, guess, sign in self.schedule_rows:
            A, c, _ = self.matrices(None, reg)
            start = state
            T = self.crossing(A, c, start, target, guess)
            end, E = self.flow(A, c, start, T)
            v = A*(end-c)
            assert not v[0, 0].contains(0), "Nontransverse crossing"
            F = E*D
            D = arb_mat([[F[i,j]-v[i,0]*F[0,j]/v[0,0] for j in range(2)] for i in range(3)])
            end[0,0] = target
            state = end
            arc = dict(region=reg, target=target, duration=T, start=start, end=end,
                       speed=v[0,0], sign=sign)
            if verify_path:
                normalized = dict(region=abs(reg), target=sign*target, duration=T,
                                  start=sign*start)
                arc["path"] = self.verify_arc(A, sign*c, normalized)
            arcs.append(arc)
        assert arcs[0]["sign"] == 1 and arcs[-1]["sign"] == -1
        zeros = [a for a in arcs if a["target"] == 0]
        assert len(zeros) == 2 and zeros[0]["speed"] < 0 and zeros[1]["speed"] > 0
        return col([state[1,0], state[2,0]]), arb_mat([[D[i,j] for j in range(2)] for i in (1,2)]), arcs


def certify(record, digits, radius="1e-35", shifted=False, parameters=None):
    f = ReturnMap(record, digits, parameters)
    center = col([arb(str(v)) for v in record["cycle"]["q"][:2]])
    for _ in range(8):
        value, J, _ = f.evaluate(center.entries())
        center = midpoint(center - (J-eye(2)).inv()*(value-center))
    if shifted:
        center[0,0] += arb("1e-12")
    box = col([arb(v, arb(radius)) for v in center.entries()])
    value0, J0, _ = f.evaluate(center.entries())
    _, J, arcs = f.evaluate(box.entries(), verify_path=True)
    C = midpoint((J0-eye(2)).inv())
    R = eye(2)-C*(J-eye(2))
    K = center-C*(value0-center)+R*(box-center)
    contraction = max(sum(abs(R[i,j]) for j in range(2)).abs_upper() for i in range(2))
    inclusion = all(box[i,0].contains_interior(K[i,0]) for i in range(2))
    if not (inclusion and contraction < 1):
        raise ValueError("Krawczyk inclusion/contraction rejected")
    trace, determinant = J.trace(), J.det()
    # Necessary and sufficient for the eigenvalues of a real 2x2 matrix to
    # lie strictly in the unit disk; applies without assuming real roots.
    schur = [1-trace+determinant, 1+trace+determinant, 1-determinant]
    assert all(s > 0 for s in schur), "Schur stability rejected"
    eig = [(trace-(trace*trace-4*determinant).sqrt())/2,
           (trace+(trace*trace-4*determinant).sqrt())/2]
    return dict(id=record["id"], passed=True, root=[bounds(v) for v in box.entries()],
                krawczyk_image=[bounds(v) for v in K.entries()], contraction=bounds(contraction),
                period=bounds(sum(a["duration"] for a in arcs)),
                multipliers=[bounds(v) for v in eig], schur_margins=[bounds(v) for v in schur],
                path_cells=[a["path"]["accepted_cells"] for a in arcs],
                flight_times=[bounds(a["duration"]) for a in arcs],
                section_itinerary=[dict(region=a["region"], target=a["target"], sign=a["sign"]) for a in arcs],
                first_upward_return=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--digits", type=int, default=80)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    source = ROOT / "results/chua-layer-atlas/four-cycles.json"
    old = json.loads(source.read_text())
    records = []
    for row in old["records"]:
        cert = certify(row, args.digits)
        records.append(cert)
        print(row["id"], cert["period"]["ball"], "certified", flush=True)
    separations = []
    for i, a in enumerate(records):
        for b in records[i+1:]:
            # Disjoint boxes for fixed points of the same FIRST return imply
            # geometrically distinct orbits, not merely distinct phases.
            ar, br = a["root"][0], b["root"][0]
            gap = abs(arb(ar["ball"])-arb(br["ball"]))
            assert gap > 0
            separations.append(dict(pair=[a["id"], b["id"]], y_gap=bounds(gap)))
    rejected = False
    try:
        certify(old["records"][0], args.digits, shifted=True)
    except ValueError as error:
        rejected = "Krawczyk" in str(error)
    assert rejected, "Shifted root control was not rejected at inclusion"
    report = dict(parameters=dict(alpha="453/50", beta="71729/5000", m0="-8/7+3/500", m1="-5/7", gamma="0"),
                  working_digits=args.digits, records=records, separations=separations,
                  shifted_root_control_rejected=rejected, passed=True,
                  scope="At least four distinct asymptotically stable periodic flow orbits. No maximum count or hiddenness claim.")
    path = args.output or OUT / f"four-cycle-certificate-{args.digits}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
