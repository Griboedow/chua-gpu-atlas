"""Recompute reported conclusions from numerical records and interval bounds.

This is a record-level cross-check, not an independent interval integrator or a
proof checker for Newhouse's theorem. It deliberately ignores stored pass flags
when recomputing signs, forecast scores, coverage and Lyapunov inequalities.
Only the requested output file is written. Use certificate --output options to
replay the underlying enclosures without overwriting the published evidence.
"""
import argparse
import copy
import hashlib
import json
import subprocess
from collections import Counter
from decimal import Decimal, localcontext
from pathlib import Path

from flint import arb, arb_mat, ctx
from chua_research_followup import ROOT, OUT, dump


def read(name):
    return json.loads((OUT / name).read_text(encoding="utf-8"))


def ball(record):
    # The short display ball can lose a strict sign through decimal rounding.
    # The union of separately enclosed endpoints remains outward rounded and
    # usually retains much more information. Never use a midpoint alone.
    if "lower" in record and "upper" in record:
        return arb(record["lower"]).union(arb(record["upper"]))
    return arb(record["ball"])


def require(condition, message):
    if not condition:
        raise ValueError(message)


def tangency(record):
    faces = record["face_signs"]
    require(ball(faces["beta_lower"]) < 0, "lower parameter face")
    require(ball(faces["beta_upper"]) > 0, "upper parameter face")
    require(ball(faces["a_lower"]) < 0, "left derivative face")
    require(ball(faces["a_upper"]) > 0, "right derivative face")
    require(ball(record["normal_curvature"]) > 0, "quadratic curvature")
    require(ball(record["parameter_unfolding"]) > 0, "parameter unfolding")
    mu, nu = map(ball, record["full_return_multipliers"])
    require(mu > 1 and nu > 0 and nu < 1, "squared-map saddle")
    det = ball(record["full_return_determinant"])
    require(det > 0 and det < 1, "squared-map dissipativity")
    published_bounds = [
        (faces["beta_lower"], "-9.334e-12", "-9.295e-12"),
        (faces["beta_upper"], "9.295e-12", "9.334e-12"),
        (faces["a_lower"], "-36.754", "-22.344"),
        (faces["a_upper"], "22.344", "36.754"),
        (record["normal_curvature"], "2.2495e16", "3.6604e16"),
        (record["parameter_unfolding"], "1.862685", "1.862993"),
        (record["full_return_multipliers"][0], "3.09870905", "3.09870925"),
        (record["full_return_multipliers"][1], "4.4316e-8", "4.4339e-8"),
        (record["full_return_determinant"], "1.3732e-7", "1.3739e-7"),
    ]
    for quantity, lower, upper in published_bounds:
        require(ball(quantity) > arb(lower) and ball(quantity) < arb(upper),
                "outward rounding in the proof tables")
    local = record["local_map"]
    D = [[ball(x) for x in row] for row in local["centered_derivative"]]
    a, b, c, d, M, Q, L, K, G = map(arb,
        ["1.7", ".000012", ".000012", ".000223", "60", ".000001",
         ".00001", "100", ".000002"])
    require(D[0][0] < -a and abs(D[0][1]) < b and abs(D[1][0]) < c
            and abs(D[1][1]) < d, "local derivative guards")
    require(ball(local["second_derivative_upper"]) < M, "state Hessian guard")
    require(ball(local["centered_parameter_derivative_bound"]) < Q,
            "centered parameter guard")
    bounds = [(c+d*L)/(a-b*L), (b+L*d)/(a-L*c), d+L*b, 1/(a-L*c),
              (M*(1+L)**3+(d+L*b)*K)/(a-b*L)**2,
              (M*(1+L)**3+K*(c*L+d)**2)/(a-L*c),
              Q*(1+L)+(d+L*b)*G, (G+Q*(1+L))/(a-L*c)]
    require(all(x < limit for x, limit in zip(bounds, [L,L,1,1,K,K,G,G])),
            "graph-transform inequalities")
    require(a-b*L > 1 and c*L+d < 1, "graph-transform domains")
    delta = arb("1e-7")
    require(abs(ball(record["initial_coordinate_interval"])) < delta,
            "departure within local domain")
    require(L*abs(ball(record["initial_coordinate_interval"])) < arb("8e-13"),
            "exact unstable graph within departure box")
    require(abs(ball(record["return_value"][1])) < delta, "arrival within graph domain")
    require(abs(ball(record["departure_witness"]["unstable_coordinate"])) > delta,
            "nontrivial excursion")
    cells = record["itinerary_cells"]
    require(len(cells) == 36 and all(len(row) == 2 and min(row) > 0 for row in cells),
            "complete return path record")
    return dict(returns=36, flights=72, accepted_time_cells=sum(map(sum,cells)),
                graph_inequalities_recomputed=10)


def thresholds(record):
    require(len(record["records"]) == 7, "seven ray records")
    L, delta = arb(".02"), arb("1e-7")
    for row in record["records"]:
        require(all(abs(ball(x)) < delta for x in row["image"]), "ray image domain")
        du, ds = map(ball, row["derivative"])
        require(du-L*abs(ds) > 0, "transverse ray")
        u,s = map(ball,row["left_image"])
        require(u+L*abs(s) < 0, "lower ray endpoint")
        u,s = map(ball,row["right_image"])
        require(u-L*abs(s) > 0, "upper ray endpoint")
        lo,hi = map(Decimal,row["bracket"])
        require(abs((hi-lo)-Decimal("2e-11")) < Decimal("1e-70"), "bracket width")
    return dict(rays=7, inequalities_per_ray=5)


def coverage(record, local):
    rows = sorted(record["records"],key=lambda r:Decimal(r["lower"]))
    lower = next(r["bracket"][0] for r in local["records"] if r["angle"] == 0)
    require(Decimal(rows[0]["lower"]) == Decimal(".299"), "coverage start")
    require(Decimal(rows[-1]["upper"]) == Decimal(lower), "coverage endpoint")
    require(all(Decimal(a["upper"]) == Decimal(b["lower"])
                for a,b in zip(rows,rows[1:])), "coverage gap or overlap")
    allowed = {"origin", "origin_by_separator", "origin_by_continuous_flow"}
    require(all(Decimal(r["lower"]) < Decimal(r["upper"]) and r["label"] in allowed
                for r in rows), "capture label or interval orientation")
    central = record["central"]
    require(Decimal(central["upper_amplitude"]) == Decimal(".299"), "central join")
    require(ball(central["maximum_absolute_x"]) < 1, "central trajectory domain")
    require(ball(central["final_energy"]) < ball(central["capture_level"]), "central capture")
    return dict(intervals=len(rows), capture_types=dict(Counter(r["label"] for r in rows)),
                limitation="Endpoint continuity and recorded destinations; replay is required to recompute every transported set.")


def neighborhoods(record):
    row, = record["records"]
    pars = row["parameters"]
    alpha,m0,m1 = (arb(str(pars[k])) for k in ["alpha","m0","m1"])
    beta = ball(row["uniform_beta_interval"])
    J = arb("29.516618475").union(arb("29.516618493"))
    require(beta.contains(J), "uniform hiddenness interval")
    require(m0 < -1 and m1 > -1, "three equilibria")
    A = arb_mat([[-alpha*(1+m1),alpha,0],[1,-1,1],[0,-beta,0]])
    entries = row["outer_lyapunov"]["P"]
    require(all(Decimal(entries[i][j]) == Decimal(entries[j][i])
                for i in range(3) for j in range(3)), "exact symmetry of Lyapunov matrix")
    P = arb_mat([[arb(x) for x in r] for r in entries])
    Q = -(A.transpose()*P+P*A)
    for matrix in [P,Q]:
        require(all(arb_mat([[matrix[i,j] for j in range(k)] for i in range(k)]).det() > 0
                    for k in [1,2,3]), "uniform Lyapunov positivity")
    xe = (m1-m0)/(1+m1)
    level = (xe-1)**2/(2*P.inv()[0,0])
    require(ball(row["captured_energy"]) < level, "uniform exit capture")
    S = arb_mat([[ball(x) for x in r] for r in row["modal_basis"]])
    W = S.inv()
    norm = sum(x*x for x in W.entries()).sqrt()
    require(arb("2.53e-7") < arb("1e-6")/(2*norm), "published origin radius")
    pnorm = max(sum(abs(P[i,j]) for j in range(3)).abs_upper() for i in range(3))
    require(arb(".00322") < (level/(2*pnorm)).sqrt(), "published outer radius")
    return dict(positive_principal_minors_recomputed=6, origin_radius="2.53e-7",
                outer_radius=".00322", limitation="Flight enclosures and the spectral root are recomputed by certificate replay.")


def predictions():
    rows=[]
    calibration_count=0
    frozen_path=OUT/"sheets/forecasts-frozen.json"
    original=subprocess.check_output(["git","show","112b01f:results/chua-followup/sheets/forecasts-frozen.json"],cwd=ROOT)
    require(original.replace(b"\r\n",b"\n") == frozen_path.read_bytes().replace(b"\r\n",b"\n"),
            "published forecast provenance")
    frozen={r["node"]:r for r in read("sheets/forecasts-frozen.json")["records"]}
    seen=set()
    for source in sorted((OUT/"sheets").glob("*-search.json")):
        data=json.loads(source.read_text())
        require(data["node"] in frozen and data["node"] not in seen, "duplicate or unknown forecast node")
        seen.add(data["node"])
        require(data["parameters"] == frozen[data["node"]]["parameters"], "search parameter identity")
        expected={r["peaks"]:r for r in frozen[data["node"]]["forecasts"]}
        require(sorted(r["peaks"] for r in data["records"]) == sorted(expected), "target family identity")
        for r in data["records"]:
            fold,sym,core = (Decimal(r[k]["coordinate"]) for k in ["fold","symmetry","core"])
            low,high = min(fold,sym),max(fold,sym)
            width=high-low
            require(width > 0 and low < core < high, "interior cycle location")
            require(abs((fold-sym)-Decimal(r["signed_width"])) < Decimal("1e-35"), "stored window width")
            tr,det = (Decimal(r["core"][k]) for k in ["trace","determinant"])
            require(1-tr+det > 0 and 1+tr+det > 0 and 1-det > 0, "numerical core stability")
            f=r["forecast"]
            require(f == expected[r["peaks"]], "forecast changed inside a search record")
            raw=low < Decimal(f["midpoint"]) < high
            corrected=low < Decimal(f["transport_midpoint"]) < high
            require(raw == r["midpoint_inside"] and corrected == r["transport_midpoint_inside"], "stored midpoint score")
            central=frozen[data["node"]]["index"] == [0,0]
            require(central == data["node"].endswith("_11"), "central calibration identity")
            if central:
                calibration_count+=1
            else:
                rows.append(dict(raw=raw,corrected=corrected,
                    fold_error=abs(Decimal(f["transport_fold"])-fold)/width,
                    width_error=abs(abs(Decimal(f["signed_width"]))/width-1)))
    require(len(rows)==48 and calibration_count==6 and seen == set(frozen), "held-out denominator")
    require(sum(r["raw"] for r in rows)==32 and sum(r["corrected"] for r in rows)==48,
            "published prediction hit counts")
    require(max(r["fold_error"] for r in rows) < Decimal(".022103") and
            max(r["width_error"] for r in rows) < Decimal(".066803"), "published prediction error bounds")
    return dict(held_out=len(rows), calibration_targets=calibration_count, raw_hits=sum(r["raw"] for r in rows),
        corrected_hits=sum(r["corrected"] for r in rows),
        maximum_corrected_fold_error_in_widths=format(max(r["fold_error"] for r in rows),".17g"),
        maximum_absolute_width_relative_error=format(max(r["width_error"] for r in rows),".17g"))


def main(output):
    ctx.dps=100
    with localcontext() as context:
        context.prec=110
        results={"predictions":predictions()}
        for digits in [140,160]:
            results[f"tangency_{digits}"]=tangency(read(f"hidden-search/hidden-tangency-{digits}.json"))
        for digits in [80,100]:
            results[f"thresholds_{digits}"]=thresholds(read(f"trigger/local-certificate-{digits}.json"))
            results[f"neighborhoods_{digits}"]=neighborhoods(read(f"hidden-search/uniform-equilibrium-exclusion-{digits}.json"))
        below=read("trigger/below-certificate-80.json")
        local=read("trigger/local-certificate-80.json")
        results["lower_ray"]=coverage(below,local)
        controls={}
        for sign in ["minus","plus"]:
            try:tangency(read(f"hidden-search/shifted-tangency-{sign}.json"))
            except ValueError as exc:controls["shifted_"+sign]=str(exc)
            else:raise ValueError("Shifted tangency control accepted")
        # Stored booleans remain true: rejection must depend on the altered data.
        corrupted=copy.deepcopy(below);del corrupted["records"][len(corrupted["records"])//2]
        try:coverage(corrupted,local)
        except ValueError as exc:controls["missing_capture_interval"]=str(exc)
        else:raise ValueError("Coverage gap accepted")
        corrupted=read("hidden-search/hidden-tangency-160.json")
        corrupted["normal_curvature"]={"ball":"[0 +/- 1]"}
        try:tangency(corrupted)
        except ValueError as exc:controls["uncertain_curvature"]=str(exc)
        else:raise ValueError("Unproved curvature accepted")
    files=list((OUT/"sheets").glob("*-search.json"))+[OUT/"sheets/forecasts-frozen.json"]
    files += [OUT/p for p in ["hidden-search/hidden-tangency-140.json","hidden-search/hidden-tangency-160.json",
        "hidden-search/uniform-equilibrium-exclusion-80.json","hidden-search/uniform-equilibrium-exclusion-100.json",
        "trigger/local-certificate-80.json","trigger/local-certificate-100.json","trigger/below-certificate-80.json",
        "hidden-search/shifted-tangency-minus.json","hidden-search/shifted-tangency-plus.json"]]
    report=dict(passed=True,results=results,negative_controls=controls,
        input_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes().replace(b"\r\n",b"\n")).hexdigest() for p in files},
        scope="Recomputed record-level scores and inequalities, with rejection controls. This does not verify the code that produced the enclosures, the analytic localization, or the imported dynamical theorem.")
    dump(output,report)
    print("Passed: 48 forecast scores, interval inequalities, coverage, and four rejection controls.")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT/".publication/review/record-check.json")
    main(parser.parse_args().output)
