"""Predict a slowly swept loss of the positive well from the frozen affine flights."""

import argparse
import json
from pathlib import Path

import mpmath as mp

from chua_fiber_switching import FrozenFiber


ROOT = Path(__file__).resolve().parents[1]
ALPHAS = (2048, 8192, 32768, 131072, 524288, 2097152)


def orbit(fiber, d):
    """Phase of the rising y=0 section and the matched x and y over a frozen turn."""
    T = fiber.period
    c, s = fiber.forcing(d, T, 1)
    R = mp.hypot(c, s)
    phase = ((mp.atan2(s, c) - mp.pi/2) % (2*mp.pi)) / fiber.omega

    def state(t):
        k, start, equilibrium = (fiber.b, 0, 0) if t <= d else (fiber.a, d, fiber.xe)
        x = (equilibrium + (1-equilibrium)*mp.exp(-k*(t-start))
             + mp.fdot(fiber.forcing_weights(k, start, t), (c, s)))
        y = c*mp.cos(fiber.omega*t) + s*mp.sin(fiber.omega*t)
        return x, y

    return R, state(phase)[0], state


def predict():
    mp.mp.dps = 50
    f = FrozenFiber(mp.mpf("0.1425"), mp.mpf(-2), -mp.mpf(5)/7)
    T = f.period
    critical = f.a*T/(f.a-f.b)
    amplitude = lambda d: f.amplitude(d, T, 1)
    r_star, x_star, state = orbit(f, critical)
    curvature = mp.diff(amplitude, critical, 2)
    product = -(f.a-f.b)**2/(2*curvature)

    def radial(t):
        x, y = state(t)
        return y*(x-y)/r_star

    drift = mp.quad(radial, [0, critical, T])
    coefficient = (drift**2/product)**(mp.mpf(1)/3)
    # Start at -8 scaled units on the stable frozen branch, fixing the
    # otherwise free Bi admixture in the Riccati solution.
    lead = mp.mpf(8)
    weight = -(mp.airyai(lead, 1) + mp.sqrt(lead)*mp.airyai(lead)) / (
        mp.airybi(lead, 1) + mp.sqrt(lead)*mp.airybi(lead))

    def zero(offset, derivative):
        return mp.findroot(lambda v: mp.airyai(-v, derivative)
                           + weight*mp.airybi(-v, derivative), offset)

    crossing, departure = zero(1, 1), zero(2.3, 0)
    rows = []
    for alpha in ALPHAS:
        scale = coefficient/mp.mpf(alpha)**(mp.mpf(2)/3)
        r0 = r_star - lead*scale
        guess = critical - mp.sqrt(2*(r0-r_star)/curvature)
        d0 = mp.findroot(lambda d: amplitude(d) - r0, guess)
        assert 0 < d0 < critical
        _, x0, _ = orbit(f, d0)
        rows.append({"alpha": alpha, "R0": float(r0), "x0": float(x0),
                     "Z0": float(f.omega*r0), "scale": float(scale),
                     "section_prediction": float(r_star+crossing*scale),
                     "two_well_prediction": float(r_star+departure*scale)})
    return {
        "parameters": {"kappa": "0.1425", "m0": "-2", "m1": "-5/7"},
        "fold": {"R": float(r_star), "x_at_y0_rising": float(x_star),
                 "R_second_derivative": float(curvature), "normal_form_product": float(product),
                 "radial_drift_per_turn_times_alpha": float(drift), "scale_prefactor": float(coefficient)},
        "initialization": "Positive stable frozen orbit at R=R_fold-8*scale, y=0 rising; initial state (x0,0,Z0) in fast variables.",
        "markers": {
            "section": "First downward passage of x through the frozen fold coordinate at rising y=0 sections. Interpolate R linearly between bracketing sections.",
            "two_well": "R at the first rising y=0 section starting a complete turn with x below -1 and above +1.",
        },
        "scaled_initial_amplitude": -8,
        "airy_section_zero": float(crossing), "airy_departure_zero": float(departure),
        "checks": {"section": "Relative excess over the fold agrees within 20% at both largest alphas.",
                   "two_well": "Relative excess agrees within 25% at the largest alpha."},
        "runs": rows,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path,
                   default=ROOT/"outputs/chua-dynamic-fold/forecast.json")
    target = p.parse_args().output
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(predict(), indent=2) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
