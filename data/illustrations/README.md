# Illustration data

Small measured datasets used by `render_animations.py`. They permit regenerating the GIFs without a GPU or the large results archive.

| File | Contents | Source in the full archive |
|---|---|---|
| `classic_trajectory.npz` | Saved time and state coordinates of the classic double scroll | `outputs/classic_trajectory.npz`, unchanged |
| `extent_candidates.npz` | All 2,249 repeated candidates; columns α, β, m₀, m₁ | `outputs/extent/verified.npz`, filtered by `(classes == 4).any(axis=1)` |
| `expanded_candidates.npz` | All 1,832 repeated adaptive candidates in the expanded domain | `outputs/expanded/verified.npz`, filtered by the same class-4 union criterion |
| `window_regular_mesh.npz` | Full vertices and triangle indices of the local near-zero/nonpositive-LLE surface | `outputs/detailed/window_regular_mesh.npz`, unchanged |

No interpolation, smoothing or random downsampling is added when preparing these files. The mesh itself is marching cubes at level 0.5 of the binary regular-node criterion, clipped at the local scan box. See the main methods document for that criterion and its limits.

The GIFs move an orthographic camera smoothly through a 60-degree horizontal arc and back over 24 seconds. Elevation stays at 23 degrees and roll at zero, so the view never tumbles or turns upside down. The sinusoidal sweep slows to a stop at each end and loops continuously. A shared palette keeps colors stable between frames. State or parameter values are never animated. The plotting box gives differently scaled coordinate ranges comparable visual lengths; Euclidean distances in the rendered image therefore are not physical distances in parameter space.

Run `python render_animations.py` from the repository root to regenerate all four GIFs. The defaults are 240 frames at 10 fps; `--frames`, `--fps` and `--sweep` control the loop duration and angular range. `docs/assets/animations.json` records the encoded duration and camera settings.
