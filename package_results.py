"""Package the full scientific results for a GitHub Release, with checksums.

Run from any directory. Extract the resulting ZIP into the repository root.
Source, illustration data and animations are already included in the Git repo.
"""
import argparse
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import hashlib
import json

ROOT=Path(__file__).resolve().parent

def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream,"sha256").hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT/"dist")
    parser.add_argument("--version",default="v1.0.0")
    args=parser.parse_args()
    if not args.version or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.-" for c in args.version):
        parser.error("Version must contain only ASCII letters, digits, dots or hyphens")
    files=[]
    for folder in [ROOT/"outputs",ROOT/"outputs"/"detailed",ROOT/"outputs"/"extent"]:
        if not (folder/"index.html").is_file():
            raise FileNotFoundError(f"Missing report: {folder.relative_to(ROOT)}/index.html")
        files.extend(p for p in folder.iterdir() if p.is_file() and p.suffix in {".npz",".json",".html",".js",".png",".pdf"})
    files.extend(ROOT/name for name in ["LICENSE","THIRD_PARTY_NOTICES.md"])
    manifest={str(p.relative_to(ROOT)).replace("\\","/"):{"bytes":p.stat().st_size,"sha256":sha256(p)} for p in sorted(files)}
    args.output.mkdir(parents=True,exist_ok=True)
    destination=args.output/f"chua-results-{args.version}.zip"
    with ZipFile(destination,"w",compression=ZIP_DEFLATED,compresslevel=6) as archive:
        for p in sorted(files): archive.write(p,p.relative_to(ROOT))
        archive.writestr("RESULTS_MANIFEST.json",json.dumps(manifest,indent=2))
    with ZipFile(destination) as archive:
        if archive.testzip() is not None: raise RuntimeError("Archive CRC verification failed")
    digest=sha256(destination)
    (args.output/"SHA256SUMS.txt").write_text(f"{digest}  {destination.name}\n",encoding="utf-8")
    print(json.dumps({"file":destination.name,"bytes":destination.stat().st_size,"sha256":digest,"files":len(files)},indent=2),flush=True)

if __name__=="__main__": main()
