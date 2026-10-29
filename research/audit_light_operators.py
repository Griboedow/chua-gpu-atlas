"""Source-equation audit for Guo et al., Nature Communications 16, 11315 (2025).

Uses the continuous heat balance (SI Eq.10), retaining the temperature at switches.
The displayed exponential responses alone omit arbitrary initial temperatures.
This audit does not reconstruct the authors' numerical implementation.
"""
import json
from pathlib import Path
import numpy as np


def single_inertial(y,theta,damping,power,threshold):
    w,v,heat=y
    return np.array([v,-damping*v+heat-w,(power*(w<threshold)-heat)/theta])


def single_relay(z,k,b,delta):
    u,v,y=z
    return np.array([v,k*(delta-u-y)-b*v,(u>0)-y])


def coupled_heat_time(state,k,b,delta_minus,delta_plus):
    """Normalized six-state continuous-heat model; away from switch surfaces.

    Each displacement and temperature is divided by its own illuminated limit.
    State x_minus,v_minus,y_minus,x_plus,v_plus,y_plus; time t/tau_heat.
    """
    xm,vm,ym,xp,vp,yp=state
    return np.array([vm,k*(ym-xm)-b*vm,(xp>delta_plus)-ym,
                     vp,k*(yp-xp)-b*vp,(xm<delta_minus)-yp])


def main():
    rng=np.random.default_rng(2026102962);errors=[]
    for _ in range(128):
        theta,damping,power,threshold=10**rng.uniform(-1,1,4)
        raw=rng.uniform(-1,2,3);w,v,heat=raw
        z=np.array([(threshold-w)/power,-theta*v/power,heat/power])
        transformed=single_inertial(raw,theta,damping,power,threshold)*np.array([-theta/power,-theta**2/power,theta/power])
        errors.append(float(np.max(abs(transformed-single_relay(z,theta**2,damping*theta,threshold/power)))))
    k=816.0170161251663;b=13.491678487182265;delta=.030904198409989477
    result=dict(source_doi='10.1038/s41467-025-66395-3',supplement='41467_2025_66395_MOESM1_ESM.pdf',
                inspected_pages=[34,35,36,37,38,39,40,41,42,43,44,45,46,47],
                equations='SI 6,10,13,17--20; continuous heat state at switching events',
                mapping='s=t/tau_heat; U=delta-w_bar/P_bar,V=-theta*w_bar_dot/P_bar,Y=w_L_bar/P_bar; k=theta^2,b=beta_bar*theta,delta=w0_bar/P_bar',
                algebra_max_abs_error=max(errors),passed=max(errors)<1e-9,
                published_single_reference=dict(theta=20,damping_bar=.8,k=400,b=16,beam_length_mm=25,threshold_mm=5,lambda_per_mW=.0026),
                existing_relay_candidate_mapping=dict(theta=float(np.sqrt(k)),damping_bar=float(b/np.sqrt(k)),delta=delta,k=k,b=b,
                    illustrative_power_mW=100,illustrative_threshold_mm=25*.0026*100*delta,
                    caveat='Mechanical/thermal ratio and damping changed from the source reference; ideal step beam, no calibrated experimental prediction'),
                novelty='Exact equivalence to an already investigated relay limit, not a third independent chaotic ODE family',
                source_ambiguity='SI exponential formulas give special step responses without arbitrary initial temperatures; a reset implementation would be a different hybrid model')
    out=Path('results/light-operators');out.mkdir(parents=True,exist_ok=True);(out/'equivalence-audit.json').write_text(json.dumps(result,indent=2))
    print(result)
    if not result['passed']:raise RuntimeError('Mapping failed')


if __name__=='__main__':main()
