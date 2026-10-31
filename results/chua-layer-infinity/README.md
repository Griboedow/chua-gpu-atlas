# Interval certificate and reproduction

This directory accompanies the [proof of infinitely many attracting windows](../../docs/chua-infinity-certificate.md). The exact Chua transect is α = 9, m₀ = −8/7, m₁ = −5/7. A quadratic homoclinic tangency is enclosed near β = 14.30635909996146.

## Recorded evidence

| File | Role |
|---|---|
| [tangency-certificate-100.json](tangency-certificate-100.json) | Main interval certificate: all 15 final gates pass |
| [tangency-certificate-80.json](tangency-certificate-80.json) | Repeat at a second working precision |
| [negative-control-plus.json](negative-control-plus.json) | β box shifted by +10⁻¹⁰; required lower-face sign fails |
| [negative-control-minus.json](negative-control-minus.json) | β box shifted by −10⁻¹⁰; required upper-face sign fails |
| [connection-center.json](connection-center.json) | Recorded numerical tangent-line root used only to locate the test box |
| [center-reconstruction.json](center-reconstruction.json) | Reconstructs that center from the earlier approximation |
| [local-certificate.json](local-certificate.json) | Initial local study and approximate basis; the main program recomputes all local bounds |
| [derivative-controls.json](derivative-controls.json) | Independent mpmath checks of local map derivatives |
| [global-derivative-controls.json](global-derivative-controls.json) | Independent 27-fold composition, curvature and unfolding checks |
| [validation-summary.json](validation-summary.json) | Expected pass/fail outcomes, unchanged forecast hash and publication checks |
| [Overview PDF](chua-infinity-overview.pdf), [interval-sign PDF](chua-infinity-signs.pdf) | Vector versions of the two report figures |

The local study alone is explicitly marked incomplete; it is not the final proof. The two negative controls are expected failures of the sufficient tangency test. All earlier forecast and orbit-check failures remain in [the prediction directory](../chua-layer-prediction/README.md).

## Run the certificate

From the repository root, use a Python environment with the CPU dependencies:

```bash
python -m pip install -r requirements-certificate.txt
python research/certify_chua_homoclinic_tangency.py --digits 80
python research/certify_chua_homoclinic_tangency.py --digits 100
```

The main run takes about four seconds on the recorded machine, excluding Python startup. The certificate uses python-flint 0.9.0 and the bundled FLINT/Arb arithmetic. Recorded versions are included in its JSON. Floating-point routines supply starting guesses only; interval inclusion and all-time path tests decide validity.

To check that the chosen parameter neighborhood matters:

```bash
python research/certify_chua_homoclinic_tangency.py --digits 100 --beta-shift 1e-10 --output results/chua-layer-infinity/negative-control-plus.json
python research/certify_chua_homoclinic_tangency.py --digits 80 --beta-shift=-1e-10 --output results/chua-layer-infinity/negative-control-minus.json
```

Both commands are expected to return a nonzero exit status after writing a failed certificate. They are not meant to prove absence of all tangencies elsewhere; they reject these particular shifted boxes.

Independent numerical controls:

```bash
python research/locate_chua_homoclinic_tangency.py
python research/check_chua_interval_derivatives.py
python research/check_chua_global_connection.py
```

The last two programs use mpmath matrix exponentials and separate crossing solves. Finite differences at two scales check local first/second derivatives and the composed normal gap, state curvature and parameter derivative. The full-composition control takes about two minutes. These controls test the implementation; they are not interval existence proofs.

To reconstruct the preliminary local study, run `python research/certify_chua_layer_organizer.py`. It regenerates the approximate coordinate basis and local enclosures. Re-run the main certificate afterwards; no unverified approximate eigenvector is treated as an exact eigendirection.

## Figures and integrity

```bash
python research/plot_chua_infinity.py
python research/manifest_refined_chua_atlas.py --check
```

The figures use the measured 15–85-peak data and the 100-digit certificate. No uncomputed windows are drawn as measurements. Error bars in the sign figure represent interval enclosures of whole faces, not statistical confidence intervals.

The manifest compares the published files byte-for-byte after documented newline normalization. Reproduction changes elapsed times and may change enclosing balls across software platforms; regenerate a manifest only when intentionally publishing a new set of evidence.

The mathematical derivation, exact scope and software trust assumptions are in the [proof report](../../docs/chua-infinity-certificate.md).
