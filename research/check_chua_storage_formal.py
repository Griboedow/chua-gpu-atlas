"""Kernel-check the storage algebra, audit axioms and reject an altered balance."""
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
from check_chua_formal import parse_axioms,ALLOWED_AXIOMS

ROOT=Path(__file__).resolve().parents[1]
PROJECT=ROOT/'formal'
SOURCE=PROJECT/'ChuaProof/Storage.lean'


def main():
    lake=shutil.which('lake')
    if not lake:raise SystemExit('lake must be on PATH; see formal/README.md.')
    source=SOURCE.read_text(encoding='utf-8')
    assert not re.search(r'\b(sorry|admit|axiom|native_decide|unsafe)\b',source)
    def run(path):
        return subprocess.run([lake,'env','lean',str(path)],cwd=PROJECT,
                              capture_output=True,text=True,encoding='utf-8',errors='replace')
    result=run(SOURCE)
    if result.returncode:raise SystemExit(result.stdout+result.stderr)
    axioms=parse_axioms(result.stdout)
    expected={'ChuaProof.'+x for x in ['storage_balance','storage_dissipation','storage_nonincrease']}
    assert set(axioms)==expected
    assert all(not (set(a)-ALLOWED_AXIOMS) for a in axioms.values())
    path=PROJECT/'.lake/build/chua-rejection-controls/WrongStorage.lean'
    path.parent.mkdir(exist_ok=True,parents=True)
    altered=source.replace('2 * p * z * phi','3 * p * z * phi',1)
    assert altered!=source
    path.write_text(altered,encoding='utf-8',newline='\n')
    rejected=run(path)
    assert rejected.returncode and 'error:' in rejected.stdout+rejected.stderr
    inputs=[SOURCE,PROJECT/'lean-toolchain',PROJECT/'lake-manifest.json',PROJECT/'lakefile.toml']
    hashes={p.relative_to(PROJECT).as_posix():hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in inputs}
    data=dict(passed=True,source_sha256_lf=hashes,theorem_axioms=axioms,
        allowed_axioms=sorted(ALLOWED_AXIOMS),altered_cross_term_rejected=True,
        axiom_output=result.stdout.strip(),negative_control_output=(rejected.stdout+rejected.stderr).strip(),
        scope='Three real-algebra implications only. Coordinate construction, Chua flow, positivity of the nonlinear integral, properness and LaSalle convergence are not formalized.')
    (ROOT/'results/chua-hidden-boundary/storage-formal-check.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('Three storage declarations kernel-checked; altered balance rejected.')


if __name__=='__main__':main()
