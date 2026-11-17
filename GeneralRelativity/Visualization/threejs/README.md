# Three.js viewer for GeneralRelativity scenes

`index.html` renders the newline-delimited JSON scene files written by
`GeneralRelativity.Visualization`. Open it directly in a browser (no server or
build step; Three.js comes from a CDN) and drop a scene file on the page, or use
the file picker. A copy of the viewer can live next to the scene files, e.g. in
`~/Documents/Visualization`:

    cp -R GeneralRelativity/Visualization/threejs/{index.html,styles} ~/Documents/Visualization/

## Visualization styles

Every scene can be shown in several styles, chosen with the **style** menu at
run time. A style that lacks the data it needs is listed as unavailable, with
the reason as a tooltip. The scene can choose the starting style and preset its
options (`Visualization.writeView`, or the `style`, `optionNames` and
`optionValues` parameters of `Visualization.SceneRecorder`).

| Style id | Shows | Needs in the scene |
|---|---|---|
| `embedding` | 3D view; the equatorial plane embedded as a surface of revolution (Flamm's paraboloid), coloured by the Kretschmann scalar; spheres; field slices | nothing (surface from `writeEmbeddingSurface`) |
| `river` | River model (Hamilton and Lisle): space flowing inward at `sqrt(1 + g_tt)`, swirling with frame dragging; light flashes expanding at c relative to the flow; a contour where the river reaches c | slice (`writeSpacetimeSlice`) |
| `river3d` | The river in 3D: concentric shells of free-falling space (latitude/longitude nets joined by radial lines) stream inward, stretch radially, cross the horizon at c and twist with frame dragging; a wedge is cut away | slice |
| `clocks` | Time-contracting grid: a sheet that sinks by the time each static clock has lost, `(1 - alpha)(t - t0)`; ticking clock faces; each body's own proper-time clock | slice |
| `spacetime` | 2+1 spacetime diagram: worldlines against coordinate time; light cones on the "now" plane that narrow and tilt with the metric | tracks; cones need the slice |
| `twins` | Twin paradox: Minkowski diagram with sliders for speed, turnaround distance and acceleration; worldlines, proper-time ticks, yearly light signals and the traveller's line of simultaneity recomputed live, with the simulated worldlines overlaid | a twin-paradox scene (`Examples.TwinParadox`) |
| `lensing` | Gravitational lensing like a deep telescope image: background galaxies seen through a foreground galaxy (isothermal ellipse with external shear) or point mass; sliders for mass, shape and redshifts; critical curves, a source-plane inset with caustics, convergence and magnification maps; drag to move the lens, scroll to zoom | a `lens` record (`Visualization.writeLens`, `Examples.GalaxyLensing`) |
| `lensfit` | Lens inversion: the observed image, the simulation of the first radius-based mass estimate, the fitted simulation and the residual in units of the noise, side by side with the fitted Einstein ring; below them the chi-squared and the mass estimate per iteration against the truth; the side panel lists every parameter before and after the fit with its standard error; scroll to change the exposure | `lensImage` records and a `lensFit` record (`Visualization.writeLensImage`, `writeLensFit`, `Examples.LensInversion`) |
| `raytrace` | GPU ray tracer of null geodesics in Kerr-Schild form: shadow, photon ring, lensed sky, thin accretion disk with Doppler and gravitational shifts, ray-traced bodies; orbit or ride along with a body | `metric` record of a known family (`writeMetric`) |

Both river styles run on free-fall time `T = t + (1/c) int (v / c_r) dr` (`v` the
river speed, `c_r` the outward radial speed of light, both from the slice); for
Schwarzschild this is Gullstrand-Painleve time. In it the river carries a body
dropped from rest far away exactly with the streaks, so streaks, shells and the
bodies (re-timed sample by sample) share one clock and true speeds; only the
light flashes of the 2D river run in slow motion. The spacetime diagram draws ct
on the same scale as length, so light far from the mass is at 45 degrees and no
worldline can look faster than its light cone; it shows a window of ct around
now that scrolls as time runs.

The river, clock and spacetime styles drop one spatial dimension: they draw the
equatorial plane; the 3D river keeps all three, using the equatorial speeds at each radius (`theta = pi/2`, or `z = 0` in a Cartesian chart). Weak fields
(Earth, cislunar) have river speeds near `4e-5 c` and clock lags near `1e-9`;
those styles exaggerate them for display and say so in the side panel.

### Showing and hiding bodies

Each body in the side panel's legend has a checkbox that hides it (with its
paths, clock and ray-traced sphere) in every style; the choice survives style
changes. A scene can start with a body hidden: `trackVisible` in
`SceneRecorder` (or `visible = false` in `writeTrack`); for example
`SchwarzschildOrbitScene(showNewtonian = false)`.

### Style options (ids for `writeView`)

| Style | Option ids (default) |
|---|---|
| `embedding` | `surface` (1), `wire` (1), `fields` (1) |
| `river` | `tracers` (2500), `flow` (1, slow-motion rate of the light flashes), `flashes` (1), `rings` (5), `map` (1) |
| `river3d` | `shells` (10), `meridians` (12), `parallels` (5), `spokes` (1), `cut` (0 none, 1 front wedge, 2 upper half), `outer` (0.7) |
| `clocks` | `warp` (1), `clocks` (1), `stride` (8), `turns` (8), `grid` (1) |
| `spacetime` | `window` (1.5, half-height of the visible ct window relative to the view), `cones` (1), `coneSize` (0.05), `rings` (6), `bodyCones` (1), `now` (1) |
| `twins` | `beta` (0.8), `distance` (4 ly), `acceleration` (1 g), `tick` (1 yr), `signals` (1), `simultaneity` (1), `simulated` (1) |
| `lensing` | `mass` (log10 solar masses), `model` (0 galaxy, 1 point), `q`, `pa`, `gamma`, `gpa`, `zl`, `zs`, `aperture` (kpc), `sources` (150, up to 400; a deep field has about 150 in the 24" sky), `stars`, `seed`, `fov` (arcsec), `lensLight`, `crit`, `ring`, `map` (0, 1 convergence, 2 magnification), `inset`, `supersample`, `exposure`; the lens position is set by dragging; the lens itself defaults to the scene's `lens` record |
| `lensfit` | `panels` (0 all, 1 observed, 2 first simulation, 3 fitted, 4 residual, 5 observed and fitted), `exposure` (1.5), `log` (1, logarithmic stretch), `range` (5, residual colour range in sigma), `ring` (1, fitted Einstein ring), `history` (1, iteration plots), `labels` (1) |
| `raytrace` | `resolution` (1, fraction of device pixels once still), `samples` (24), `steps` (600), `accuracy` (0.04), `disk` (1), `diskOuter` (8, in rs), `shift` (1), `sky` (0 stars, 1 grid, 2 checker), `bodySize` (0.35, in rs), `camera` (0 orbit, k ride with body k), `paths` (0), `exposure` (1) |

Booleans are 0 or 1. For example, a scene that opens in the ray tracer riding
with its first body:

    GeneralRelativity.Visualization.SceneRecorder scene(..., style = "raytrace",
      optionNames = {"camera", "sky"}, optionValues = {1, 2});

### Ray-tracer refinement

While the view or the bodies move, the ray tracer draws a coarse preview (30% of
the window). Once everything is still for 0.2 s it traces a sharper preview,
then full device resolution with finer integration steps, then keeps adding
sub-pixel-jittered samples (up to `samples`) that average into an anti-aliased
image. Full-resolution passes are traced in horizontal bands sized to a
per-frame pixel budget, so the page stays responsive; the side panel shows the
progress. The integrator uses an analytic gradient of the Kerr-Schild
Hamiltonian (one metric evaluation per Runge-Kutta stage), and exposure is
applied at display time, so changing it does not re-trace.

### Ray-tracer metrics

The ray tracer integrates Hamilton's equations for `H = g^{mu nu} p_mu p_nu / 2`
with `g = eta + f k k` (Kerr-Schild form, regular across the horizon; the
inverse is exactly `eta - f K K`), following the approach of
<https://michaelmoroz.github.io/TracingGeodesics/>. Metric names it knows:
`minkowski`, `schwarzschild` (`p = {rs}`), `kerr` (`{rs, a}`),
`reissnerNordstrom` (`{rs, rQ}`), `kerrNewman` and `kerrSchild` (`{rs, a, rQ}`).
The same formula is `GeneralRelativity.Metrics.kerrSchild`, checked against the
Boyer-Lindquist metrics in `Tests.GeneralRelativityVisualizationTests`. The
image is a snapshot: bodies are drawn where they are now, without light-travel
delay, and the camera is a static observer.

## Adding a style

A style is one plain script in `styles/` (not an ES module, so it loads from
disk) that appends a definition to `window.GRStyles`; list it with a `<script>`
tag in `index.html`. The definition has `id`, `label`, `order`, `description`,
`options` (sliders, checkboxes and selects that the core turns into controls),
`available(data)` (null, or the reason it cannot run) and `create(ctx)`, which
returns an object with:

- `group`: a `THREE.Group` the core adds to the scene;
- `place(x, sample, t)`: where a body at scene coordinates `x` is drawn; the core
  draws bodies, trails and the legend with it (`dynamicPlacement: true` re-places
  them every frame; `hidesBodies: true` hides the core's bodies);
- optional `hidesTracks` (true, or a getter, to hide the simulated tracks), `update(clock, dt)`, `render(renderer, camera, scene)` (to draw the
  frame itself, as the ray tracer does), `camera()`, `readout()`,
  `onOption(id, value)` (return true when applied in place) and `dispose()`.

`ctx` provides the parsed scene (`data`: header, slice and its fields, tracks,
spheres, metric, embedding rows, profile), the clock, a colour map and colour
bar, a bilinear `sampler(slice, field, mask)` for slice fields, coordinate
helpers, and `setTimeRange(t0, t1)` for a style that computes its own motion
(the twin paradox does, so its sliders can change the trip's length).

## Scene records

`scene` (`beginScene`), `view` (`writeView`), `metric` (`writeMetric`),
`surfaceRow` (`writeEmbeddingSurface`), `profile` (`writeRadialProfile`), `field`
(`writeFieldSlice`), `slice` and `sliceField` (`writeSpacetimeSlice`,
`writeSliceField`), `sphere` (`writeSphere`), `track` (`writeTrack`) and
`sample` (`writeSample`), `lens` (`writeLens`), `lensImage` and `lensFit`
(`writeLensImage`, `writeLensFit`). `Visualization.SceneRecorder` writes the
trajectory records from one component. Samples written without a `track` belong
to one default track, so older scene files still load.

Produce files with `GeneralRelativity.Examples.SchwarzschildOrbitScene`
(`schwarzschild_scene.ndjson`: the Schwarzschild geodesic and the Newton-Cartan
geodesic), `KerrOrbitScene` (`kerr_scene.ndjson`), `TwinParadox` (`twin_paradox.ndjson`), `GalaxyLensing` (`galaxy_lensing.ndjson`),
`LensInversion` (`lens_inversion.ndjson`, plus `lens_inversion_observed.mat`, the synthetic
observation that `LensInversionFromFile` re-analyses into `lens_inversion_file.ndjson`;
`tools/lens_image_to_mat.py` makes such a file from a real picture), `EarthJ2OrbitScene`
(`earth_scene.ndjson`) or `CislunarOrbitScene` (`cislunar_scene.ndjson`). They
are written to `~/Documents/Visualization` when `~/Documents` exists; set the
`sceneDirectory` parameter or the `GR_SCENE_DIR` environment variable to use
another folder, otherwise they go to the simulation working directory.
