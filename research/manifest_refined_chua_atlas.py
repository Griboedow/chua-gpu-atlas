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
continue_chua_slope_families check_chua_slope_families map_chua_slope_coexistence plot_chua_slope_families
resolve_chua_overlap_geometry check_chua_overlap_geometry plot_chua_overlap_geometry
predict_chua_layers chua_precision_flights check_chua_precision_flights verify_chua_layer_predictions
audit_chua_layer_predictions describe_chua_layer_words analyze_chua_prediction_evidence
trace_chua_layer_organizer plot_chua_layer_predictions
chua_interval_map chua_interval_connection certify_chua_layer_organizer
certify_chua_homoclinic_tangency locate_chua_homoclinic_tangency
check_chua_interval_derivatives check_chua_global_connection plot_chua_infinity
audit_chua_coexistence check_chua_formal plot_chua_proof_audit
certify_chua_coexisting_cycles certify_chua_hidden_cycle study_chua_stability
check_chua_accessibility prepare_chua_accessibility probe_chua_accessibility
check_chua_stability_prediction study_chua_circuit_constraints plot_chua_accessibility
refine_chua_slope_atlases check_chua_slope_refinement study_chua_slope_grazing
check_chua_slope_branches plot_chua_slope_refinement
derive_chua_global_sector study_chua_hidden_boundary certify_chua_basin_saddle
check_chua_hidden_boundary check_chua_storage_formal plot_chua_hidden_boundary
certify_chua_basin_separator check_chua_basin_proof audit_chua_local_results plot_chua_local_results'''.split()


def digest(path):
    raw=path.read_bytes()
    if path.suffix in ['.py','.cl','.txt','.md','.json','.html','.csv','.js','.cjs','.svg','.lean','.toml','.yml','.ps1'] or path.name == 'lean-toolchain':
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
    paths.update(p for p in (ROOT/'results/chua-overlap-geometry').iterdir() if p.is_file())
    paths.update(p for p in (ROOT/'results/chua-layer-prediction').iterdir() if p.is_file())
    paths.update(p for p in (ROOT/'results/chua-layer-infinity').iterdir() if p.is_file())
    paths.update(p for p in (ROOT/'results/chua-proof-audit').iterdir() if p.is_file())
    paths.update(p for p in (ROOT/'results/chua-accessibility').iterdir() if p.is_file())
    paths.update(p for p in (ROOT/'results/chua-slope-refinement').iterdir() if p.is_file())
    paths.update(p for p in (ROOT/'results/chua-hidden-boundary').iterdir() if p.is_file())
    paths.update(p for p in (ROOT/'results/chua-basin-proof').iterdir() if p.is_file())
    paths.update(ROOT/p for p in ['docs/chua-basin-and-coexistence.md',
        'docs/assets/chua-basin-proof.png','docs/assets/chua-coexistence-localization.png'])
    paths.add(ROOT/'docs/chua-hidden-boundary.md')
    paths.update(ROOT/'docs/assets'/name for name in ['chua-global-sector.png',
        'chua-basin-saddle.png','chua-linear-contact.png'])
    paths.add(ROOT/'docs/chua-slope-refinement.md')
    paths.update(ROOT/'docs/assets'/name for name in ['chua-slope-refined-volumes.png',
        'chua-slope-refined-sections.png','chua-slope-grazing-law.png'])
    paths.update(ROOT/p for p in ['docs/chua-accessibility.md','docs/assets/chua-accessibility.png',
        'docs/assets/chua-four-cycle-certificate.png','docs/assets/chua-circuit-accessibility.png'])
    paths.update(p for p in (ROOT/'formal').glob('*') if p.is_file() and p.name != '.gitignore')
    paths.update((ROOT/'formal/ChuaProof').glob('*.lean'))
    paths.update(ROOT/p for p in ['docs/chua-coexistence-audit.md','docs/assets/chua-proof-audit.png','research/check_chua_formal_local.ps1'])
    paths.update(ROOT/p for p in ['requirements-certificate.txt','docs/chua-infinity-certificate.md','docs/research-guide.md','docs/README.md','research/README.md'])
    paths.update((ROOT/'docs/assets').glob('chua-infinity-*.png'))
    paths.add(ROOT/'requirements-precision.txt')
    paths.add(ROOT/'docs/chua-layer-prediction.md')
    paths.update((ROOT/'docs/assets').glob('chua-layer-prediction-*.png'))
    paths.add(ROOT/'docs/assets/chua-layer-motion-description.png')
    paths.add(ROOT/'docs/assets/chua-layer-predicted-motions.png')
    paths.add(ROOT/'research/check_chua_slope_volume_viewer.cjs')
    paths.add(ROOT/'research/check_chua_overlap_guide.cjs')
    paths.update([ROOT/'README.md',ROOT/'docs/chua-refined-atlas.md',ROOT/'docs/chua-layer-census.md'])
    paths.add(ROOT/'docs/chua-family-gallery.md')
    paths.add(ROOT/'docs/chua-slope-plane.md')
    paths.add(ROOT/'docs/chua-slope-volumes.md')
    paths.add(ROOT/'docs/chua-slope-families.md')
    paths.add(ROOT/'docs/chua-overlap-geometry.md')
    paths.add(ROOT/'results/README.md')
    paths.update((ROOT/'docs/assets').glob('chua-overlap-*.png'))
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
        text_hash_normalization='CRLF is converted to LF for .py, .cl, .txt, .md, .json, .html, .csv, .js, .cjs, .svg, .lean, .toml, .yml, .ps1 and lean-toolchain; binary files are unchanged.',
        scope='Hashes identify the finite catalogue, forecasts, interval certificate and partial formal proof. Integrity checks do not themselves prove mathematical claims; see docs/chua-infinity-certificate.md and docs/chua-coexistence-audit.md for the arguments and remaining obligations.')
    manifest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('Hashed',len(paths),'files')


if __name__=='__main__':main()
