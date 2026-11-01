"""Build the Lean pilot, audit transitive axioms, and exercise two rejection controls."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
THEOREMS = ["not_min_at_left", "not_min_at_right", "tangency_box", "quadratic_contact", "graph_guards"]
ALLOWED_AXIOMS = {"propext", "Classical.choice", "Quot.sound"}


def digest(path):
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def parse_axioms(output):
    found = {}
    for name, body in re.findall(r"'([^']+)' depends on axioms:\s*\[([^\]]*)\]", output):
        found[name] = sorted(x.strip() for x in body.split(",") if x.strip())
    for name in re.findall(r"'([^']+)' does not depend on any axioms", output):
        found[name] = []
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, default=ROOT / "formal",
                        help="Optional identical build copy, useful for a different drive.")
    parser.add_argument("--output", type=Path, default=ROOT / "results/chua-proof-audit/formal-check.json")
    args = parser.parse_args()
    project = args.project.resolve()
    lake = shutil.which("lake")
    if not lake:
        raise SystemExit("lake is not on PATH; install the pinned Lean toolchain first.")
    inputs = [Path("lean-toolchain"), Path("lakefile.toml"), Path("lake-manifest.json"),
              Path("ChuaProof.lean"), Path("Audit.lean")]
    inputs += sorted(p.relative_to(ROOT/"formal") for p in (ROOT/"formal/ChuaProof").rglob("*.lean"))
    hashes = {}
    for relative in inputs:
        original = ROOT/"formal"/relative
        hashes[relative.as_posix()] = digest(original)
        if digest(project/relative) != digest(original):
            raise SystemExit(f"Build copy differs: {relative}")
        if relative.suffix == ".lean" and re.search(
                r"\b(sorry|admit|axiom|native_decide|unsafe)\b", original.read_text(encoding="utf-8")):
            raise SystemExit(f"Disallowed proof shortcut in {relative}")

    def run(*command, success=True):
        proc = subprocess.run([lake, *command], cwd=project, text=True,
                              encoding="utf-8", errors="replace", capture_output=True)
        if success and proc.returncode:
            raise SystemExit(proc.stdout + proc.stderr)
        return proc

    build = run("build")
    print("Lean build passed", flush=True)
    audit = run("env", "lean", "Audit.lean")
    dependencies = parse_axioms(audit.stdout)
    expected = {"ChuaProof." + name for name in THEOREMS}
    if set(dependencies) != expected:
        raise SystemExit("Unexpected axiom report: " + audit.stdout)
    if any(set(values) - ALLOWED_AXIOMS for values in dependencies.values()):
        raise SystemExit("Nonstandard dependency in the proof: " + audit.stdout)

    controls = project/".lake/build/chua-rejection-controls"
    controls.mkdir(parents=True, exist_ok=True)
    altered = (project/"ChuaProof/TangencyBox.lean").read_text(encoding="utf-8")
    altered = altered.replace("0 < g a bR", "g a bR < 0")
    reversed_top = controls/"ReversedTop.lean"
    reversed_top.write_text(altered, encoding="utf-8")
    bad = run("env", "lean", str(reversed_top), success=False)
    rejected_sign = bad.returncode != 0 and "type mismatch" in (bad.stdout + bad.stderr).lower()
    if not rejected_sign:
        raise SystemExit("Reversed upper-face control did not fail as expected: " + bad.stdout + bad.stderr)
    extra = controls/"ExtraAxiom.lean"
    extra.write_text("axiom unsupportedClaim : False\n"
                     "theorem unsafeConclusion : False := unsupportedClaim\n"
                     "#print axioms unsafeConclusion\n", encoding="utf-8")
    extra_output = run("env", "lean", str(extra))
    extra_deps = parse_axioms(extra_output.stdout)
    rejected_axiom = bool(set(extra_deps.get("unsafeConclusion", [])) - ALLOWED_AXIOMS)
    if not rejected_axiom:
        raise SystemExit("The extra-axiom control was not detected.")
    version = run("env", "lean", "--version").stdout.strip()
    result = dict(
        status="passed", toolchain=version,
        source_sha256_lf=hashes, theorem_axioms=dependencies,
        allowed_axioms=sorted(ALLOWED_AXIOMS),
        controls=dict(reversed_parameter_face_rejected=rejected_sign,
                      extra_axiom_rejected=rejected_axiom),
        build_output=build.stdout.strip(), axiom_output=audit.stdout.strip(),
        scope="Generic real-analysis tangency-box implication and exact rational guards. "
              "No formal Chua flow, interval evaluator, invariant-manifold construction, "
              "or Newhouse--Robinson theorem is included.",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print("All five declarations use only standard axioms; both rejection controls passed.", flush=True)


if __name__ == "__main__":
    main()
