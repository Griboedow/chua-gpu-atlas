"""Hash or verify the published refinement, its code, and prior evidence inputs."""
import argparse
import ast
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]
SCRIPTS='''audit_refined_chua_atlas bridge_chua_itineraries build_chua_family_registry
check_chua_event_map check_fast_flights check_sheet_interpolation chua_basin_capture
chua_event_map chua_fast_flights chua_overlap_basins chua_window_ancestry
export_refined_chua_atlas extend_chua_adding manifest_refined_chua_atlas plot_refined_chua_atlas plot_chua_family_gallery
refine_chua_adding_check refine_chua_satellites refine_chua_sheets scan_chua_offcenter
select_chua_offcenter verify_chua_offcenter verify_four_chua_cycles
refine_readme_layers refine_readme_thin_cores plot_readme_layers
scan_chua_slope_plane check_chua_slope_plane plot_chua_slope_plane
scan_chua_slope_volumes check_chua_slope_volumes plot_chua_slope_volumes build_chua_slope_volume_viewer
continue_chua_slope_families check_chua_slope_families map_chua_slope_coexistence plot_chua_slope_families'''.split()


def digest(path):
    raw=path.read_bytes()
    if path.suffix in ['.py','.cl','.txt','.md','.json','.html','.csv','.js','.cjs']:
        raw=raw.replace(b'\r\n',b'\n')
    return hashlib.sha256(raw).hexdigest()


def code_dependencies():
    pending=[ROOT/'research'/(name+'.py') for name in SCRIPTS];found=set()
    while pending:
        path=pending.pop()
        if path in found:continue
        found.add(path)
        for node in ast.walk(ast.parse(path.read_text(encoding='utf-8'))):
            names=[a.name for a in node.names] if isinstance(node,ast.Import) else [node.module] if isinstance(node,ast.ImportFrom) and node.module else []
            for name in names:
                for folder in [ROOT,ROOT/'research']:
                    local=folder/(name.split('.')[0]+'.py')
                    if local.exists():pending.append(local)
    return found | {ROOT/'kernels.cl',ROOT/'requirements.txt'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args();out=ROOT/'results/chua-layer-atlas';manifest=out/'manifest.json'
    if args.check:
        data=json.loads(manifest.read_text());bad=[]
        for path,expected in data['sha256'].items():
            if not (ROOT/path).is_file() or digest(ROOT/path)!=expected:bad.append(path)
        if bad:raise SystemExit('Checksum mismatch: '+', '.join(bad))
        print('Verified',len(data['sha256']),'published, code, and input hashes');return
    paths=code_dependencies() | {p for p in out.iterdir() if p.is_file() and p!=manifest}
    paths.update(p for p in (out/'readme-layer-view').rglob('*') if p.is_file())
    paths.update(p for p in (ROOT/'results/chua-slope-plane').iterdir() if p.is_file())
    paths.update(p for p in (ROOT/'results/chua-slope-volumes').iterdir() if p.is_file())
    paths.update(p for p in (ROOT/'results/chua-slope-families').iterdir() if p.is_file())
    paths.add(ROOT/'research/check_chua_slope_volume_viewer.cjs')
    paths.update([ROOT/'README.md',ROOT/'docs/chua-refined-atlas.md',ROOT/'docs/chua-layer-census.md'])
    paths.add(ROOT/'docs/chua-family-gallery.md')
    paths.add(ROOT/'docs/chua-slope-plane.md')
    paths.add(ROOT/'docs/chua-slope-volumes.md')
    paths.add(ROOT/'docs/chua-slope-families.md')
    paths.update((ROOT/'docs/assets/chua-families').glob('*.png'))
    paths.update(ROOT/'docs/assets'/name for name in ['chua-refined-atlas.png','chua-refined-adding.png','chua-four-cycle-overlap.png'])
    paths.update(ROOT/'docs/assets'/name for name in ['chua-readme-layers.gif','chua-readme-layers-poster.png','chua-readme-layer-zoom.png'])
    paths.update(ROOT/'docs/assets'/name for name in ['chua-slope-plane.png','chua-slope-peak-patterns.png'])
    paths.update(ROOT/'docs/assets'/name for name in ['chua-slope-volumes.png','chua-slope-volumes-interior.png','chua-slope-volume-sections.png','chua-slope-volumes.gif'])
    paths.update(ROOT/'docs/assets'/name for name in ['chua-slope-family-sheets.png','chua-slope-family-sheets.gif','chua-slope-coexistence.png','chua-slope-basin-sensitivity.png','chua-slope-itinerary-change.png'])
    paths.update(ROOT/'results/chua-layer-census'/name for name in ['family-edges.json','census-analysis.json','cycles.json','cycles-extra.json',
        'adding-cycles.json','second-surfaces.json','second-boundaries.json','atlas-comparison.npz'])
    paths.add(ROOT/'results/chua-layers/extended.json')
    data=dict(created_utc=datetime.now(timezone.utc).isoformat(),base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        sha256={p.relative_to(ROOT).as_posix():digest(p) for p in sorted(paths)},
        text_hash_normalization='CRLF is converted to LF for .py, .cl, .txt, .md, .json, .html, .csv, .js and .cjs; binary files are unchanged.',
        scope='Hashes identify this finite numerical catalogue and its executable evidence, not mathematical certification.')
    manifest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('Hashed',len(paths),'files')


if __name__=='__main__':main()
