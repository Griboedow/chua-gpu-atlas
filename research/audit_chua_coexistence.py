"""Repeat the tangency gates and record their scope in a coexistence audit.

This program checks numerical premises. It does not implement or prove the
Newhouse--Robinson theorem or the localization of its generated sinks.
"""

import hashlib
import json
from pathlib import Path

from certify_chua_homoclinic_tangency import certify

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/chua-proof-audit"


def main():
    OUT.mkdir(exist_ok=True)
    runs = []
    for digits, shift in [(80, "0"), (100, "0"), (100, "1e-10"), (80, "-1e-10")]:
        data = certify(digits, shift)
        expected = "passed" if shift == "0" else "failed"
        assert data["status"] == expected, (digits, shift, data["status"])
        runs.append(dict(
            digits=digits, beta_shift=shift, status=data["status"],
            expected_status=expected, gates=data["gates"],
            local_gates=data["local_gates"], graph_gates=data["graph_gates"],
            full_return_multipliers=data["full_return_multipliers"],
            full_return_determinant=data["full_return_determinant"],
            face_signs=data["face_signs"], normal_curvature=data["normal_curvature"],
            parameter_unfolding=data["parameter_unfolding"],
            global_flights=sum(len(row) for row in data["itinerary_cells"]),
            global_time_cells=sum(sum(row) for row in data["itinerary_cells"]),
            software=data["software"],
        ))
        print(f"digits={digits}, shift={shift}: {data['status']} (expected {expected})", flush=True)
    source_files = [
        "research/certify_chua_homoclinic_tangency.py",
        "research/certify_chua_saddle.py",
        "research/chua_interval_map.py", "research/chua_interval_connection.py",
        "results/chua-layer-infinity/connection-center.json",
        "results/chua-layer-infinity/local-certificate.json",
    ]
    evidence = dict(
        status="local_numerical_hypotheses_reproduced",
        coexistence_status="conditional_on_localization_argument",
        formal_scope="Generic tangency-box lemma and rational guards only; see formal/.",
        runs=runs,
        source_sha256={p: hashlib.sha256((ROOT/p).read_bytes().replace(b"\r\n", b"\n")).hexdigest()
                       for p in source_files},
        open_obligation=("Supply a precise local Newhouse--Robinson theorem or a parameter-dependent "
                         "compact extension, with the generated sinks confined to the unchanged "
                         "Chua return-map neighborhood. This audit does not mark that step as proved."),
    )
    (OUT/"numerical-audit.json").write_text(json.dumps(evidence, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
