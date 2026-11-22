"""Forecast the first exit from a saddle neighborhood after voltage precharge."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from observed_saddle_flights import NOMINAL, Section


OUT = Path(__file__).resolve().parents[1] / "results/chua-saddle-residence"
BREAK_VOLTAGE = 1.1606
SECTION_SEED = [-2.321, -0.787, 3.129]
RADIUS = 0.02
OFFSETS = (0.008, 0.004, 0.002, 0.001, 0.0003, 0.0001, 0.00003,
           0.00001, 0.000003, 0.000001, 0.0000003, 0.0000001)


def save(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def locate(section, point, left, depth):
    def displacement(x):
        _, state, derivative = section.launch(x)
        yz, dx = state[1:], derivative[1:, 0]
        for _ in range(depth):
            _, yz, dp = section.full_return(yz)
            dx = dp @ dx
        return float(left @ (yz-point)), float(np.linalg.norm(yz-point)), float(left @ dx)
    a = brentq(lambda x: displacement(x)[0], 7.14, 7.16, xtol=4e-15)
    gap, distance, slope = displacement(a)
    return dict(depth=depth, x=a, projected_displacement=gap,
                section_distance=distance, projected_launch_derivative=slope)


def make_forecast():
    section = Section()
    point, period, eigenvalues, left = section.saddle()
    lam = eigenvalues[-1]
    crossings = [locate(section, point, left, depth) for depth in range(4)]
    center = crossings[2]["x"]
    # Divide out the two full unstable returns to obtain sensitivity at entry.
    gain = crossings[2]["projected_launch_derivative"] / lam**2
    predictions = []
    for offset in OFFSETS:
        for side in (-1, 1):
            fractional = np.log(RADIUS / (offset * gain)) / np.log(lam)
            predictions.append(dict(offset_x=side*offset, launch_x=center+side*offset,
                                    predicted_fractional_return=float(fractional),
                                    predicted_first_exit=int(np.ceil(fractional))))
    record = dict(source="https://doi.org/10.1007/s11071-022-08078-y",
                  parameters=NOMINAL, published_seed=SECTION_SEED,
                  section="x=-1 with increasing x; one crossing per full saddle period",
                  point=point.tolist(), period=float(period), multipliers=list(eigenvalues),
                  left_unstable_unit=left.tolist(), crossings=crossings, crossing_x=center,
                  breakpoint_V=BREAK_VOLTAGE, crossing_V=center*BREAK_VOLTAGE,
                  entry_gain_per_x=gain, projected_exit_radius=RADIUS,
                  count_definition="First section index with absolute projected displacement >= radius; first hit has index zero.",
                  interpolation="Linear in log absolute displacement between the last two section returns.",
                  prediction="log(radius/(abs(offset_x)*entry_gain_per_x))/log(unstable_multiplier); round upward for the integer first exit.",
                  tests={"main":"All 24 integer forecasts within one full return, and maximum fractional error below 0.35 return.",
                         "independent":"Full unsplit DOP853 at the eight offsets of magnitude 0.008, 0.001, 0.0001 and 0.00001: same integer counts.",
                         "negative_control":"Replace unstable multiplier with 1.5 at unchanged gain and exit radius; compare fractional errors.",
                         "slope_control":"The increase from 0.001 to 0.0001 and from 0.0001 to 0.00001 should differ from log(10)/log(lambda) by <0.20 turn on each side."},
                  predictions=predictions)
    save("forecast.json", record)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["forecast"])
    args = parser.parse_args()
    if args.stage == "forecast":
        make_forecast()
