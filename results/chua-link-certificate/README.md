# Certified S02/S08 link

The [research note](../../docs/chua-link-certificate.md) gives the parameters, isotopy argument and Jones polynomial obstruction. The [standard run](certificate.json) and [finer mesh](refinement.json) record the dyadic coordinate grid, vertex counts and digests, coordinate error bounds in grid units, box comparisons and exact diagram polynomials. Polygon digests use consecutive little-endian signed 64-bit integer triples.

Run from the repository root with the certificate, precision and knot dependencies. On versions of Python supported by the pinned Numba and Spherogram packages, these commands write to `outputs/chua-link-certificate/`.

```bash
python -m pip install -r requirements-certificate.txt -r requirements-precision.txt -r requirements-knots.txt
python research/certify_chua_knot.py --figure
python research/certify_chua_knot.py --digits 160 --steps-per-time 160 --output outputs/chua-link-certificate/refinement.json
```

Passing `--steps-per-time 16` exercises the failure control. It rejects S02 segments 24 and 376 on a 466-edge polygon. During each successful run the verifier reconstructs the two original 140-digit section boxes, matches every serialized field and repeats the Krawczyk inclusion. The selected geometry precision applies to the subsequent orbit enclosures. Each projected crossing is interior and transverse, with a strictly ordered pair of depths and a unique position along its edge. The joint polynomials match the [earlier exhaustive bracket calculation](../chua-link-search/overlaps/S02-S08/exhaustive-jones.json). [Verifier](../../research/certify_chua_knot.py)
