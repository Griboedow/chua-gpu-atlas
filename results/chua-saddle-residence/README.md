# Residence near the Chua saddle cycle

The [study report](../../docs/chua-saddle-residence.md) describes the voltage-only launch and the full-return exit convention. The model uses the mathematical column of Table 2 in [Kuznetsov et al.](https://doi.org/10.1007/s11071-022-08078-y).

| File | Contents |
|---|---|
| [forecast.json](forecast.json) | Published before measuring the 24 target launches. Includes the orbit, projected launch derivative, voltage threshold, offset pairs, predicted counts and acceptance criteria. |
| [scores.csv](scores.csv) | One row per starting voltage with predicted and observed exits. |
| [measurements.json](measurements.json) | Full section crossings, two event-bracketing steps, independent DOP853 controls, two post-exit examples, and scores. |
| [voltage-start-residence.png](voltage-start-residence.png) | Voltage launch, saddle returns and forecast errors. |

From the repository root, run `python research/observed_saddle_validate.py --figure`. The regenerated JSON, CSV and PNG go to `outputs/chua-saddle-residence/`. The forecast and its original generator are preserved in [commit 1de6255](https://github.com/Griboedow/chua-gpu-atlas/commit/1de6255e4dd35a9dbf8af5a8564d853b90b6dd64).
