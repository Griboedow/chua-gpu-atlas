"""Twelve targeted modal transverse controls, with per-case CPU time bounds.

Reuses the literal absolute-value drag and mode_control from the structural
audit. Never overwrites that earlier audit. No GPU or full-pair trajectories.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import math
from pathlib import Path
import subprocess
import sys
import time
import warnings

import numpy as np

from audit_rayleigh_pair import mode_control

ROOT = Path(__file__).resolve().parents[1]


def sanitize(value):
    if isinstance(value, dict):
        return {k: sanitize(v) for k, v in value.items()}
    if isinstance(value, list):
        return [sanitize(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def annotate_small_multipliers(result):
    """Keep raw matrix eigenvalues but suppress exponents of unresolved modes."""
    period = result["period_sigma"]
    for data in result["integrators"].values():
        phi = np.asarray(data["endpoint"])[6:10].reshape(2, 2)
        raw = np.asarray(data["transverse_multipliers"])
        moduli = np.hypot(raw[:, 0], raw[:, 1])
        threshold = max(1e-12, 1e-12*np.linalg.norm(phi, ord=2))
        unresolved = np.flatnonzero(moduli <= threshold).tolist()
        data["unresolved_transverse_multiplier_indices"] = unresolved
        data["multiplier_resolution_heuristic_threshold"] = float(threshold)
        data["raw_transverse_exponents_from_matrix"] = data["transverse_exponents"].copy()
        for i in unresolved:
            data["transverse_exponents"][i] = None
        rho = float(max(moduli))
        if unresolved and rho > 100*threshold:
            log_leading = float(np.log(rho))
            log_other = float(-data["endpoint"][-1]-log_leading)
            data["trace_inferred_log_abs_multipliers"] = [log_leading, log_other]
            data["trace_inferred_transverse_exponents"] = [log_leading/period, log_other/period]
            data["trace_inference_scope"] = "Small log modulus inferred from integrated trace and the resolved dominant multiplier; not a direct FP64 eigenvalue measurement."
    return result


def worker(mu, kappa, symmetric):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = mode_control(mu, kappa, symmetric)
    result["warnings"] = [str(w.message) for w in caught]
    for method, data in result["integrators"].items():
        endpoint = np.asarray(data["endpoint"])
        phi = endpoint[6:10].reshape(2, 2)
        norm = np.linalg.norm(phi, ord=2)
        logdet = float(-endpoint[-1])
        expected = float(np.exp(logdet))
        resolvable = expected > 1e-12 * max(norm*norm, np.finfo(float).tiny)
        data["liouville_log_abs_determinant"] = logdet
        data["direct_determinant_relative_check_resolvable"] = bool(resolvable)
        data["direct_determinant_relative_error"] = float(abs(np.linalg.det(phi)-expected)/expected) if resolvable else None
        data["small_multiplier_warning"] = "Do not interpret zero/tiny rounded multipliers or a small absolute determinant residual as relative accuracy. Use the integrated trace as the determinant reference."
        data["state_closure"] = float(max(abs(endpoint[0]), abs(endpoint[1]-result["section_velocity"])))
    d, r = result["integrators"]["DOP853"], result["integrators"]["Radau"]
    reliable = (result["dop853_radau_full_augmented_difference"] < 1e-6
                and max(d["flow_multiplier_error"], r["flow_multiplier_error"]) < 1e-6
                and max(d["state_closure"], r["state_closure"]) < 1e-7)
    lo = min(d["maximum_transverse_multiplier"], r["maximum_transverse_multiplier"])
    hi = max(d["maximum_transverse_multiplier"], r["maximum_transverse_multiplier"])
    result["status"] = ("transversely unstable numerical modal cycle" if reliable and lo > 1+1e-6
                        else "transversely stable numerical modal cycle" if reliable and hi < 1-1e-6
                        else "unresolved at the stated checks")
    result["scope"] = "A modal transverse instability, if present, would not itself establish chaos."
    return sanitize(annotate_small_multipliers(result))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--mu", type=float)
    parser.add_argument("--kappa", type=float)
    parser.add_argument("--symmetric", action="store_true")
    parser.add_argument("--wall-limit", type=float, default=420)
    parser.add_argument("--case-limit", type=float, default=60)
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(args.mu, args.kappa, args.symmetric), allow_nan=False))
        return
    started = time.perf_counter()
    cases = [(10., kappa, symmetric) for kappa in (.01, 1.5, 4.) for symmetric in (True, False)]
    cases += [(30., kappa, symmetric) for kappa in (.1, 1.5, 4.) for symmetric in (True, False)]
    output = {"audit_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
              "equations": "x_i''=-x_i-kappa*(x1+x2)+mu*(1-abs(x_i'))*x_i'",
              "scope": "Twelve hand-chosen modal controls only, no generic 4D trajectories or parameter grid.",
              "time_unit": "canonical sigma", "physical_scope": "Strong-mu points are mathematical controls, not calibrated LCE operating conditions.",
              "script": "research/rayleigh_modal_followup.py",
              "base_solver": "research/audit_rayleigh_pair.py:mode_control",
              "requested_cases": [{"mu": m, "kappa": k, "mode": "symmetric" if s else "antisymmetric"} for m, k, s in cases],
              "results": [], "wall_limit_seconds": args.wall_limit, "case_limit_seconds": args.case_limit}
    target = ROOT/"results/rayleigh-pair/modal-followup.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    for mu, kappa, symmetric in cases:
        remaining = args.wall_limit-(time.perf_counter()-started)
        if remaining < 1:
            output["termination"] = "wall limit reached; incomplete cases were not started"
            break
        command = [sys.executable, "-X", "utf8", str(Path(__file__).resolve()), "--worker", "--mu", str(mu), "--kappa", str(kappa)]
        if symmetric:
            command.append("--symmetric")
        case_start = time.perf_counter()
        try:
            process = subprocess.run(command, capture_output=True, text=True, encoding="utf8", timeout=min(args.case_limit, remaining))
            if process.returncode:
                entry = {"mu": mu, "kappa": kappa, "mode": "symmetric" if symmetric else "antisymmetric", "status": "numerical worker failed", "error": process.stderr[-3000:]}
            else:
                entry = json.loads(process.stdout)
        except subprocess.TimeoutExpired:
            entry = {"mu": mu, "kappa": kappa, "mode": "symmetric" if symmetric else "antisymmetric", "status": "CPU time limit; unresolved"}
        entry["wall_seconds"] = time.perf_counter()-case_start
        output["results"].append(entry)
        output["wall_seconds"] = time.perf_counter()-started
        target.write_text(json.dumps(output, indent=2, allow_nan=False)+'\n', encoding='utf8')
        print({"mu": mu, "kappa": kappa, "mode": entry["mode"], "status": entry["status"],
               "rho": entry.get("integrators", {}).get("DOP853", {}).get("maximum_transverse_multiplier"),
               "case_seconds": entry["wall_seconds"]}, flush=True)
    else:
        output["termination"] = "all twelve controls completed"
    output["wall_seconds"] = time.perf_counter()-started
    target.write_text(json.dumps(output, indent=2, allow_nan=False)+'\n', encoding='utf8')
    print('Saved', target, 'wall', output['wall_seconds'], flush=True)


if __name__ == '__main__':
    main()
