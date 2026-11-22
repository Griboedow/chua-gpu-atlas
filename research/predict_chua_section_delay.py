"""Compute the first finite-alpha correction to the fold-coordinate passage marker."""

import argparse
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from chua_fiber_switching import FrozenFiber
from predict_chua_fold_delay import orbit


ROOT = Path(__file__).resolve().parents[1]
ALPHAS = (16384, 65536, 262144, 1048576, 4194304)


def feedback(fiber, d, r_star, x_star):
    """Compute finite-epsilon return derivative E and radial derivative D_x by affine flows."""
    w, a, b, T = map(float, (fiber.omega, fiber.a, fiber.b, fiber.period))
    c, s = map(float, fiber.forcing(d, fiber.period, 1))
    t0 = (np.arctan2(s,c)-np.pi/2) % (2*np.pi) / w
    # Fixed phase is the ascending y=0 crossing within the inner flight.
    assert 0 < t0 < float(d)
    # Variables are x,y,Z,1; their epsilon variations X,Y,W; then y*p,Z*p,
    # and the accumulated integral of y*p/R. Here p'= -g'(x)*p, p(0)=1.
    v = np.array([x_star, 0, w*r_star, 1, 0, 0, 0, 0, w*r_star, 0.], dtype=float)
    for k, const, dt in [(b, 0, float(d)-t0), (a, float(fiber.xe)*a, T-float(d)), (b, 0, t0)]:
        M = np.zeros((10,10))
        M[0,:4] = [-k, 1, 0, const]
        M[1,2], M[2,1] = 1, -w*w
        M[4,4:6] = [-k, 1]
        M[5,[0,1,6]] = [1,-1,1]
        M[6,5] = -w*w
        M[7,7:9], M[8,7:9] = [-k,1], [-w*w,-k]
        M[9,7] = 1/r_star
        v = expm(M*dt) @ v
    assert max(abs(v[:3]-[x_star,0,w*r_star])) < 1e-11
    # At the section, y'=wR and x'=-b*x. Correct the perturbed return time.
    E = v[4] - (-b*x_star)*v[5]/(w*r_star)
    return float(E), float(v[9])


def coefficients():
    mp.mp.dps = 50
    f = FrozenFiber(mp.mpf(".1425"),mp.mpf(-2),-mp.mpf(5)/7)
    T = f.period
    d = f.a*T/(f.a-f.b)
    amp = lambda v: f.amplitude(v,T,1)
    xsec = lambda v: orbit(f,v)[1]
    r_star, x_star, state = orbit(f,d)
    r2,r3 = mp.diff(amp,d,2),mp.diff(amp,d,3)
    x1,x2 = mp.diff(xsec,d),mp.diff(xsec,d,2)
    q1 = f.a-f.b
    B = q1/(2*x1)
    A = -2*B*x1*x1/r2
    C = (-A*r3-4*B*x1*x2-x1*q1*q1)/(2*x1*r2)
    D = (q1*q1-2*B*x2-C*r2)/(6*x1*x1)

    def radial(t):
        x,y = state(t)
        return y*(x-y)/r_star

    drift = float(mp.quad(radial,[0,d,T]))
    E, drift_x = feedback(f,d,float(r_star),float(x_star))
    A,B,C,D = map(float,(A,B,C,D))
    c,d3,e,h = C/(A*B), D/B**2, E/(A*drift), drift_x/(B*drift)

    def field(s, u):
        z, correction = u
        dz = s+z*z
        ddz = 1+2*z*dz
        return [dz, 2*z*correction + c*z*s + d3*z**3 + e-h*z*dz-ddz/2]

    def marker(_, u):
        return u[0]

    s0 = -8.
    z0 = -np.sqrt(-s0)
    correction0 = -(c*z0*s0 + d3*z0**3)/(2*z0)
    sol = solve_ivp(field,(s0,1.2),[z0,correction0],method="DOP853",
                    rtol=2e-13,atol=2e-14,events=marker,max_step=.04)
    if not sol.success or len(sol.t_events[0]) != 1:
        raise ArithmeticError("Riccati correction marker did not converge")
    zero = sol.t_events[0][0]
    shift = -sol.y_events[0][0,1]/zero
    return f, {"R_fold":float(r_star),"x_fold":float(x_star),"critical_inner_duration":float(d),
               "A":A,"B":B,"C":C,"D":D,"E":E,"radial_drift":drift,
               "radial_drift_x":drift_x,"Airy_marker":zero,
               "leading_coefficient":(drift**2/(A*B))**(1/3)*zero,
               "correction_coefficient":drift*shift,
               "initial_value_U1":correction0, "marker_value_U1":float(sol.y_events[0][0,1])}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,default=ROOT/"outputs/chua-dynamic-fold/section-forecast.json")
    target=p.parse_args().output
    f,k = coefficients()
    rows=[]
    critical=mp.mpf(k["critical_inner_duration"])
    curvature=mp.diff(lambda v:f.amplitude(v,f.period,1),critical,2)
    for alpha in ALPHAS:
        leading=k["leading_coefficient"]*alpha**(-2/3)
        scale=leading/k["Airy_marker"]
        r0=k["R_fold"]-8*scale
        guess=critical-mp.sqrt(-16*scale/curvature if curvature>0 else 16*scale/(-curvature))
        d0=mp.findroot(lambda v:f.amplitude(v,f.period,1)-r0,guess)
        _,x0,_=orbit(f,d0)
        rows.append({"alpha":alpha,"R0":r0,"x0":float(x0),"Z0":float(f.omega)*r0,
                     "leading_excess":leading,"corrected_excess":leading+k["correction_coefficient"]/alpha})
    record={"status":"Follow-up to the failed leading-order forecast. No measured values enter its coefficients.",
            "markers":"Only the interpolated x section marker from forecast.json is scored.",
            "coefficients":k,"runs":rows,
            "check":"At the three largest unseen alphas the corrected error must be below 10% of the positive Airy leading excess and smaller than the leading-only error."}
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8",newline="\n")


if __name__=="__main__":
    main()
