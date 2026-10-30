# Physical meminductor: follow-up sources and an inaccessible thesis

Source audit, 30 October 2026. No new simulations.

**This follow-up did not recover a calibrated magnetic or optical closure for the Dinavahi–Harris moving-coil experiment.** It did identify an earlier primary source not recorded in our previous audit: Dinavahi's 2023 doctoral dissertation. Its repository metadata are accessible through the search index, but its full text was not accessible during this check. The dissertation remains a specific outstanding source, rather than evidence for or against a chaotic attractor.

## What was checked

The starting point was the existing [energy audit](meminductor-energy-audit.md), [optical-feedback report](photofeedback-research.md) and [supplement record](../results/meminductor-energy/supplement-audit.json). The official supplement had already been read; neither its retrieval nor the 3,072-trajectory illustrative survey was repeated.

| Follow-up | Actual access | Result |
| --- | --- | --- |
| Crossref record for the 2025 article [1] | Live metadata API | The returned record contains a preprint relation and no `update-to` entry. This limited metadata result does not certify the absence of corrections. |
| OpenAlex forward citations for `W4416389366` | Live work record and complete returned two-item citing list | The same MOSFET emulator paper and hardware-security review previously checked; no new indexed citing work in this snapshot. |
| Dinavahi dissertation [2] | Indexed official repository record and abstract; ordinary HTTP and Edge visits attempted | New bibliographic lead. Full text unavailable: the repository returned a Cloudflare verification page with HTTP 403. No challenge bypass was attempted. |
| Mondal and Chandra, 2026 [3] | Publisher abstract and references | A CMOS emulator, not a moving coil. Its reference list cites the 2023 physical-device paper. It does not resolve the present apparatus's closure. The paywalled full text was not examined. |
| Spiking-activity emulator paper, 2026 [4] | Publisher-indexed abstract; direct page retrieval timed out | Reports chaotic and periodic spiking in a circuit combining memristive and meminductive emulators. This is a different physical implementation; full equations were not inspected. |

Targeted searches combined the article DOI, author names, correction/erratum terms, permanent-magnet flux, coenergy, optical feedback, and the dissertation title. Search coverage and failed accesses are retained in the [machine-readable record](../results/candidate-audits/meminductor-followup-audit.json). These are bounded checks, not an exhaustive priority review.

## Why the dissertation matters, and what remains unknown

The Texas A&M record dates the dissertation to 10 October 2023 and lists `DINAVAHI-DISSERTATION-2023.pdf`, approximately 5.89 MB. Its abstract describes physical realization, small-signal local activity, passive coupling and persistent oscillations. It predates the 2025 paper and could contain useful derivations or experimental details. We did not obtain its equations, appendices, component calibration or raw measurements. Therefore we cannot say whether it supplies the missing permanent-magnet flux linkage, reciprocal force/flux law, LED/LDR transfer, response time or friction law.

In the accessible 2025 article, *edge of chaos* names a local small-signal condition: a stable biased element can invert fluctuations at some frequencies. The reported optical-circuit examples are periodic oscillations; more complicated dynamics are discussed as possible future extensions. This terminology is not a Lyapunov-spectrum demonstration of a bounded strange attractor. Conversely, the absence of such a demonstration in the inspected material does not establish that the apparatus can never exhibit chaos.

The decision remains **physical-model closure unresolved**. No additional state, delay or magnetic law has been fitted or invented here. The conditional fixed-DC energy exclusion and the negative illustrative optical survey retain their original, separate scopes. A useful next step requires ordinary access to the dissertation or a published calibrated constitutive model; another arbitrary parameter scan would not answer the experimental question.

## References

1. **Dinavahi, A., and Harris, H. R.** (2025). Edge-of-chaos operation and persistent dynamics for neuromorphic meminductor computing. *Scientific Reports* **15**, 40672. [DOI](https://doi.org/10.1038/s41598-025-12529-y). [Open text](https://pmc.ncbi.nlm.nih.gov/articles/PMC12630749/).
2. **Dinavahi, A.** (2023). *Meminductor—Physical Realization and Edge of Chaos Operation*. Doctoral dissertation, Texas A&M University. [Official repository record](https://oaktrust.library.tamu.edu/items/4dc0c47e-a888-4e98-8d78-0aacecae55f3). [Persistent handle](https://hdl.handle.net/1969.1/202865). Metadata and abstract inspected; PDF not accessed.
3. **Mondal, K., and Chandra, S.** (2026). High-frequency flux-controlled meminductor emulator with non-volatile memory. *Microsystem Technologies* **32**, 160. [DOI](https://doi.org/10.1007/s00542-026-06131-z). Publisher abstract and reference list inspected.
4. **Spiking activity in a meminductive and memristive emulator-based bionic circuit** (2026). *Chinese Physics B*. [Publisher record](https://cpb.cpsjournals.org/EN/10.1088/1674-1056/ae4c72). [DOI](https://doi.org/10.1088/1674-1056/ae4c72). Publisher-indexed abstract only; accepted manuscript online 3 March 2026 according to that record.
