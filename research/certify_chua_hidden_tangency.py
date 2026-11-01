"""Apply the shared exact-invariant-graph certificate to the hidden-lobe return."""
import argparse
import json
from chua_research_followup import OUT,dump
from chua_interval_lobe import IntervalLobe,saddle_center
from certify_chua_homoclinic_tangency import certify


def main(digits,shift,ra,rb,rv):
    source=OUT/'hidden-search/connection-center-14-2.json';center=json.loads(source.read_text())
    config=dict(source=source,map_factory=lambda d:IntervalLobe(center['parameters'],d),saddle_center=saddle_center,
        point_radius='1e-100',inflation_floor='1e-115',ra=ra,rb=rb,rv=rv,
        centered_local_derivative=True,moving_frame=IntervalLobe.moving_frame,
        guards=dict(a='1.7',b='.000012',c='.000012',d='.000223',M='60',P='.000001',L='.00001',K='100',Gamma='.000002'),
        parameter_sign='1',curvature_sign='1',
        map_description='Ordinary positive-lobe first return to x=1 with positive x velocity; two flights (outer, central). Squaring gives positive saddle multipliers.')
    try:
        result=certify(digits,shift,config)
        dump(OUT/'hidden-search'/f'hidden-tangency-{digits}.json',result)
        print(json.dumps({k:result[k] for k in ['status','gates','normal_curvature','parameter_unfolding','face_signs','return_value']},indent=2),flush=True)
    except (ValueError,AssertionError) as exc:
        dump(OUT/'hidden-search'/f'hidden-tangency-{digits}-failure.json',dict(error=str(exc),configuration={k:v for k,v in config.items() if k not in ['map_factory','saddle_center','source']}))
        raise


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--digits',type=int,default=140);ap.add_argument('--shift',default='0');ap.add_argument('--ra',default='1e-15');ap.add_argument('--rb',default='5e-12');ap.add_argument('--rv',default='8e-13');a=ap.parse_args();main(a.digits,a.shift,a.ra,a.rb,a.rv)
