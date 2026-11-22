# Capacitor precharge records

The [report](../../docs/chua-precharge.md) gives the model, measured component ranges, launch protocol and numerical criteria. The 99 settings use a five-segment diode. All launches start with zero inductor current.

| File | Fields and units |
|---|---|
| [settings.csv](settings.csv) | Row index joins to the outcome columns. R, R1–R6 and RL are in ohms, C1 and C2 in farads, L in henries, B_volts in volts. The three model coefficients, diode slopes and outer breakpoint are dimensionless. Blank components denote parameter-only settings. |
| [outcomes.csv](outcomes.csv) | Each row is one pair of initial capacitor voltages in volts. Each case column gives the finite-time outcome for the matching settings index. Signs indicate the equilibrium or mean inner-oscillation coordinate. |
| [survey.json](survey.json) | Nominal voltage bands, grid points with agreement across settings and transverse voltages, integration step and horizon in dimensionless time, and counts for step refinement. |
| [sensitivity.json](sensitivity.json) | Each pair brackets a numerical change of outcome in volts. Coupled and table-box rows follow settings order. Component sensitivities are in mV per +0.1% relative change. |
| [controls.json](controls.json) | Rectangle grid outcomes, adaptive integrations and comparison of diode models. Tail x and section spread are dimensionless. Multiply x by the breakpoint voltage B and period by R C2 for volts and seconds. |
| [transfer.json](transfer.json) | Model coefficients, threshold forecasts in volts and errors in mV. Separate integrations refine all 64 test thresholds. DOP853 records verify the sides of three predicted thresholds. |
| [voltage-preparation.png](voltage-preparation.png) | Computed outcome bands and predictions for the coupled component design. |

## Reproduce on a CPU

The recorded environment used Python 3.12.13, NumPy 2.5.5, SciPy 1.18.1, Numba 0.68.1 and Matplotlib 3.11.2. Numba compiles the RK4 batch solver on first use. From the repository root, run all four commands with the same output directory. The comparison reads the survey and sensitivity records there.

```bash
python research/chua_precharge.py --part survey --output-dir outputs/chua-precharge
python research/chua_precharge.py --part sensitivity --output-dir outputs/chua-precharge
python research/chua_precharge.py --part controls --output-dir outputs/chua-precharge
python research/chua_precharge.py --part comparison --output-dir outputs/chua-precharge
```
