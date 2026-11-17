# Three.js viewer for GeneralRelativity scenes

`index.html` renders the newline-delimited JSON scene files written by
`GeneralRelativity.Visualization`. Open it directly in a browser (no server or
build step; Three.js comes from a CDN) and drop a scene file on the page, or use
the file picker. A copy of the viewer can live next to the scene files, e.g. in
`~/Documents/Visualization`:

    cp -R GeneralRelativity/Visualization/threejs/{index.html,render.js,audio.js,render_server.py,styles} ~/Documents/Visualization/

## Visualization styles

Every scene can be shown in several styles, chosen with the **style** menu at
run time. A style that lacks the data it needs is listed as unavailable, with
the reason as a tooltip. The scene can choose the starting style and preset its
options (`Visualization.writeView`, or the `style`, `optionNames` and
`optionValues` parameters of `Visualization.SceneRecorder`).

| Style id | Shows | Needs in the scene |
|---|---|---|
| `embedding` | 3D view; the equatorial plane embedded as a surface of revolution (Flamm's paraboloid), coloured by the Kretschmann scalar; spheres; field slices. A merging binary: the orbital plane sunk into the two masses' wells and rippled by the outgoing gravitational waves (two-armed spiral at the retarded time) | nothing (surface from `writeEmbeddingSurface`; a binary: its `metric` record and the `h+`, `hx` scalars) |
| `river` | River model (Hamilton and Lisle): space flowing inward at `sqrt(1 + g_tt)`, swirling with frame dragging; light flashes expanding at c relative to the flow; a contour where the river reaches c | slice (`writeSpacetimeSlice`) |
| `river3d` | The river in 3D: concentric shells of free-falling space (latitude/longitude nets joined by radial lines) stream inward, stretch radially, cross the horizon at c and twist with frame dragging; a wedge is cut away | slice |
| `clocks` | Time-contracting grid: a sheet that sinks by the time each static clock has lost, `(1 - alpha)(t - t0)`; ticking clock faces; each body's own proper-time clock | slice |
| `spacetime` | 2+1 spacetime diagram: worldlines against coordinate time; light cones on the "now" plane that narrow and tilt with the metric | tracks; cones need the slice |
| `twins` | Twin paradox: Minkowski diagram with sliders for speed, turnaround distance and acceleration; worldlines, proper-time ticks, yearly light signals and the traveller's line of simultaneity recomputed live, with the simulated worldlines overlaid | a twin-paradox scene (`Examples.TwinParadox`) |
| `lensing` | Gravitational lensing like a deep telescope image: background galaxies seen through a foreground galaxy (isothermal ellipse with external shear) or point mass; sliders for mass, shape and redshifts; critical curves, a source-plane inset with caustics, convergence and magnification maps; drag to move the lens, scroll to zoom | a `lens` record (`Visualization.writeLens`, `Examples.GalaxyLensing`) |
| `lensfit` | Lens inversion, one case or a library: a gallery of cases with thumbnails and how close each fit came to the truth; for the chosen case its key numbers (mass with error, ring-radius estimate, truth, Einstein radius, fit quality), the observed image, the simulation of the first radius-based estimate, the fitted simulation and the residual in units of the noise, each with a colour bar and the fitted critical curve; the cost and mass per iteration; a table of every parameter before and after the fit with its error and the truth. Hover an image for pixel values; arrow keys step through cases. Drawn as HTML over the view | `lensImage` records and `lensFit` records, tagged with a case (`Visualization.writeLensInversion`; `Examples.LensInversion`, `LensInversionLibrary`, `LensInversionFromFile`) |
| `raytrace` | GPU ray tracer of null geodesics in Kerr-Schild form: shadow, photon ring, lensed sky, thin accretion disk with Doppler and gravitational shifts, ray-traced bodies; orbit or ride along with a body | `metric` record of a known family (`writeMetric`) |
| `cosmos` | The expanding universe in 3D: comoving galaxies carried apart by a(t), coloured by whether their light has reached us and whether our light can still reach them; the Hubble sphere, the particle and event horizons and the last-scattering sphere; an inset with a(t) and the energy budget; proper or comoving coordinates; log-time clock | an `expansion` record (`Examples.CosmicExpansion`) |
| `horizons` | Spacetime diagram of the expanding universe (after Davis and Lineweaver 2004): our past light cone (the teardrop), the Hubble sphere, particle and event horizons, comoving galaxies; proper or comoving distance against cosmic or conformal time | an `expansion` record |
| `cmb` | The cosmic microwave background: the model's temperature power spectrum against the Planck 2018 data and best fit, with residuals; a random sky with the model spectrum (Mollweide map or sphere, with an l range); the sound waves of the plasma in real space (the expanding acoustic shell) and in k; recombination | a `cmbSpectrum` record (`Examples.CosmicMicrowaveBackground`) |

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
| `embedding` | `surface` (1), `wire` (1), `fields` (1); a binary: `waves` (1), `waveGain` (1.5, the merger wave rises this times 6% of the region), `falloff` (0 shows `r h`, 1 the true `1/r`), `region` (100 M), `depth` (2.5) |
| `river` | `tracers` (2500), `flow` (1, slow-motion rate of the light flashes), `flashes` (1), `rings` (5), `map` (1) |
| `river3d` | `shells` (10), `meridians` (12), `parallels` (5), `spokes` (1), `cut` (0 none, 1 front wedge, 2 upper half), `outer` (0.7) |
| `clocks` | `warp` (1), `clocks` (1), `stride` (8), `turns` (8), `grid` (1) |
| `spacetime` | `window` (1.5, half-height of the visible ct window relative to the view), `cones` (1), `coneSize` (0.05), `rings` (6), `bodyCones` (1), `now` (1) |
| `twins` | `beta` (0.8), `distance` (4 ly), `acceleration` (1 g), `tick` (1 yr), `signals` (1), `simultaneity` (1), `simulated` (1) |
| `lensing` | `mass` (log10 solar masses), `model` (0 galaxy, 1 point), `q`, `pa`, `gamma`, `gpa`, `zl`, `zs`, `aperture` (kpc), `sources` (150, up to 400; a deep field has about 150 in the 24" sky), `stars`, `seed`, `fov` (arcsec), `lensLight`, `crit`, `ring`, `map` (0, 1 convergence, 2 magnification), `inset`, `supersample`, `exposure`; the lens position is set by dragging; the lens itself defaults to the scene's `lens` record |
| `lensfit` | `case` (0, index into the scene's cases), `panels` (0 all four, 1 observed and fitted, 2 observed, fitted and residual), `exposure` (1.2), `log` (1, logarithmic stretch), `range` (4, residual colour range in sigma), `ring` (1, fitted critical curve and centres), `plots` (1, iteration plots), `table` (1, parameter table) |
| `raytrace` | `resolution` (1, fraction of device pixels once still), `samples` (24), `steps` (600), `accuracy` (0.04), `disk` (1 around a mass, 0 in flat spacetime), `diskOuter` (8, in rs), `shift` (1), `sky` (0 stars, 1 grid, 2 checker), `bodySize` (0.25, in rs), `camera` (0 orbit, k ride with body k), `paths` (0), `exposure` (1) |
| `cosmos` | `frame` (0 proper, 1 comoving), `logTime` (1), `galaxies` (6000), `spheres` (1), `inset` (1) |
| `horizons` | `panel` (0 proper distance vs time, 1 comoving vs time, 2 comoving vs conformal time), `tMax` (30 Gyr), `dMax` (60 Gly), `galaxies` (1) |
| `cmb` | `panel` (0 spectrum, 1 Mollweide sky, 2 sphere, 3 sound waves, 4 recombination), `data` (1), `theory` (1), `components` (0), `skyL` (512), `lMin` (2), `lMax` (1024), `seed` (7), `rMax` (220 Mpc) |

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
progress.

With the camera still, only the bodies change, so they are traced as a
separate layer. The refined image is traced without bodies and keeps refining
while the animation plays. A mask pass (half resolution) marks the pixels whose
rays pass through a tube around the stretch of path each body is on: a 3D
occupancy grid around the recorded track from just behind the body to a few rs
ahead, widened by the body radius and the integration step. Each frame
re-traces only those pixels, with the bodies, and lays the result over the
refined image (one sample at an adaptive resolution while they move,
accumulating up to `samples` once they stop). The mask is redrawn when a body
nears the end of its stretch, and the stretch is lengthened or shortened so
that this happens about every second. Moving the camera, riding with a body or
changing an option re-traces everything as before. The integrator uses an analytic gradient of the Kerr-Schild
Hamiltonian (one metric evaluation per Runge-Kutta stage), and exposure is
applied at display time, so changing it does not re-trace.

### Lensing performance

The lensed image is drawn by one fragment shader that ray-shoots each pixel
back to the source plane through the lens equation. The galaxies sit in a
64 x 64 spatial index over the source plane (each cell lists the galaxies
whose cut-off radius reaches it), so a pixel evaluates only the few galaxies
near where its ray lands, and the cost hardly depends on `sources`. While the
lens is dragged or a slider moves, the image is drawn with one sample per
pixel; 150 ms after the last change it is redrawn once with 2 x 2
supersampling. Caustics are computed for the lens at the origin and translated
with it, so dragging does not recompute them.

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

### A merging binary in 3D

For a `binary` scene the embedding style draws the orbital plane as a sheet,
recomputed on the GPU every frame:

- **Wells:** the sheet sinks by the clocks' lag `1 - alpha` of the superposed
  metric of the two masses where they are now (the remnant after the merger).
- **Horizons and stars:** black inside a horizon, white inside a star.
- **Waves:** the outgoing gravitational wave adds a ripple, orange where
  positive and blue where negative:

      h(r, phi, t) = (M/r) 1/2 Re[(h+ - i hx)(t - r/c) e^{2 i phi}]

  This is the leading-order wave in the orbital plane, from the face-on
  waveform the scene records (the tracks' `h+` and `hx` scalars, `R h / M`),
  at the retarded time. It forms a two-armed spiral moving out at c, chirping
  as the orbit shrinks, then the merger and ringdown burst. Before the first
  sample the waveform's phase is continued backwards, so no false wave front
  leaves the origin at the start.
- **Display scaling:** the near zone, inside about 1.6 times the largest
  orbital radius, is faded out, because the formula holds only far from the
  source. The ripple is exaggerated (`waveGain`) and shown as `r h` unless
  `falloff = 1`.

The wavelength is about 160 M early in the inspiral and about 18 M at the
merger, so the sheet extends well beyond the orbit: `region`, 100 M by
default.

### Binary metric family

`binary` (`Examples.BinaryMerger`, `p = {m1, m2, R1, R2, tMerge, Mf, aF, Rremnant}`)
is a compact binary whose two bodies are the scene's first two tracks. Before
`tMerge` (in `ct`) the metric is the superposed Kerr-Schild metric of the two
masses at their current positions, `g = eta + sum (2 m_a / r_a) k_a k_a` with
`k_a = (1, n_a)` (`Metrics.binaryKerrSchild`; an approximation, not a solution
of the field equations). Its inverse has no closed form, so the shader inverts
the 4 x 4 metric and uses `dp_i/dlambda = u^mu u^nu d_i g_mu nu / 2` with
`u = g^-1 p`. `R_a > 0` makes a body a star (a neutron star) instead of a black
hole (horizon at `r_a = 2 m_a`). From `tMerge` on, the remnant is a Kerr hole of
mass `Mf` and spin length `aF` about `+z` (a star of radius `Rremnant` if that
is positive). The scene writes this record twice, and the viewer uses the
last: once at the start with `tMerge = 1e30`, and again at the merger with the
actual time. Every frame of a binary is re-traced in full, because the holes
lens everything.

## The sound of a merger

For a scene whose tracks carry `h+` and `hx` scalars (`Examples.BinaryMerger`),
the side panel adds **Gravitational-wave sound**:

- **Sound profile:** the face-on strain `h+` and the gravitational-wave
  frequency against time to the merger. A playhead follows the clock, and a
  second one follows the audio while it plays.
- **Play chirp:** the true sound. Its time scale is the scene's `gwSource`
  record: seconds per `M` for the binary's mass in solar masses. The default
  black holes sweep from 17 Hz to a 275 Hz ringdown in 0.59 s; the neutron
  stars from 290 Hz to the kHz range in 75 ms.
- **Slower versions** (4×, 16×): keep the pitch. The sound is resynthesized
  from the strain's amplitude and unwrapped phase, rather than played back
  faster or slower.
- **Pitch ×2 to ×8:** for speakers that cannot play the lowest notes.
- **Save .wav:** writes the chosen version as 16-bit mono at 44.1 kHz.
- **Live tone:** a tone that follows the animation, at the wave's frequency
  and loudness at the clock.

## Cosmology

`Examples.CosmicExpansion` (`cosmic_expansion.ndjson`, opens in `cosmos`) and
`Examples.CosmicMicrowaveBackground` (`cmb.ndjson`, opens in `cmb`) write these
records:

- `cosmology`: parameters and derived numbers;
- `expansion`: tables against cosmic time, in Gly, Gyr and c = 1;
- `recombination`, `cmbSpectrum`, `cmbData`, `cmbTheory` and `cmbComparison`;
- `acoustic`: transfer functions on a k grid at snapshot times.

How the CMB views use them:

- **Data:** the CMB style draws the Planck 2018 data (Planck Collaboration
  2020, A&A 641, A1 and A6, from the Planck Legacy Archive) next to the model.
- **Sky map:** synthesized in the browser with orthonormal spherical
  harmonics. The recurrences run per ring, with an FFT around each ring; up to
  `l = 512` it takes about 0.2 s.
- **Sound-wave view:** radially Fourier transforms the transfer functions:
  `delta(r) = (1/2 pi^2) int k^2 T(k) sin(kr)/(kr) dk`. That gives the
  profile around a point-like initial perturbation (Eisenstein, Seo and White
  2007).
- **Clock:** the CMB style runs its own clock over the snapshots
  (`ownClock`).

## Rendering a video

**Render…** in the top bar renders the animation frame by frame:

- **Camera:** stays where it is.
- **Clock:** steps through a time range (by default the whole run) at a frame
  rate (24 fps by default) and a video length that you choose.
- **Waiting:** each frame waits until the active style reports it final. For
  the ray tracer, that means every refinement sample (**samples per frame**)
  plus the bodies' layer.
- **Resolution:** set in the dialog, independent of the window. The view shows
  the frame letterboxed while it renders.

Outputs:

| Format | Made by | Saved |
|---|---|---|
| MP4 (H.264) or WebM (VP9) | the browser (WebCodecs, with `mp4-muxer` / `webm-muxer` from jsDelivr) | through the browser's download |
| PNG frames in a `.zip` | the browser | through the browser's download |
| MPEG-2 `.mpeg`, H.264 `.mp4` (ffmpeg), or PNG frames | the helper `render_server.py` | in the output folder typed in the dialog |

With a browser download, set the browser to ask where to save to choose the
folder:

- Firefox: Settings > Downloads > "Always ask you where to save files".
- Safari: Settings > General > File download location > "Ask for each
  download".
- Chrome and Edge: offer a file picker (**Choose…**).

To write straight into a folder, and to get `.mpeg` files, run the helper and
open the viewer from it:

    cd ~/Documents/Visualization && python3 render_server.py
    # then open http://127.0.0.1:8765/ and load a scene as usual

The helper:

- serves the viewer folder;
- writes the frames into the output folder;
- runs `ffmpeg`, which must be on the PATH, then deletes the frames unless
  **also keep the PNG frames** is ticked.

MPEG-2 allows only 23.976, 24, 25, 29.97, 30, 50, 59.94 and 60 fps.

Safety:

- The helper listens on 127.0.0.1 only and tells the browser to revalidate
  every file, so an updated viewer is never run from a stale cache.
- It accepts writes only from pages it served, using a per-run token, and
  checks the Host header.
- It writes only inside your home folder.

## Adding a style

A style is one plain script in `styles/` (not an ES module, so it loads from
disk) that appends a definition to `window.GRStyles`; list it with a `<script>`
tag in `index.html`. The definition has `id`, `label`, `order`, `description`,
`options` (sliders, checkboxes and selects that the core turns into controls),
`available(data)` (null, or the reason it cannot run) and `create(ctx)`, which
returns an object with:

- `group`: a `THREE.Group` the core adds to the scene;
- optional `ownClock: true` for a style that animates without tracks (it sets
  its range with `ctx.setTimeRange`);
- `place(x, sample, t)`: where a body at scene coordinates `x` is drawn; the core
  draws bodies, trails and the legend with it (`dynamicPlacement: true` re-places
  them every frame; `hidesBodies: true` hides the core's bodies);
- optional:
  - `hidesTracks` (true, or a getter, to hide the simulated tracks);
  - `update(clock, dt)`;
  - `render(renderer, camera, scene)` (to draw the frame itself, as the ray
    tracer does);
  - `camera()` and `readout()`;
  - `onOption(id, value)` (return true when applied in place);
  - `converged()` (a progressive style: true once the current frame is final;
    the Render dialog waits for it; without it a frame is taken after two
    draws);
  - `dispose()`.

`ctx` provides the parsed scene (`data`: header, slice and its fields, tracks,
spheres, metric, embedding rows, profile), the clock, a colour map and colour
bar, a bilinear `sampler(slice, field, mask)` for slice fields, coordinate
helpers, and `setTimeRange(t0, t1)` for a style that computes its own motion
(the twin paradox does, so its sliders can change the trip's length).
2D plots: `styles/canvas2d.js` (`GRPlot`) draws with the Canvas 2D API on a
canvas the size of the drawing buffer. It shows the result as a full-screen
texture, so plots are part of the WebGL image and of offline renders. It
provides axes, lines and legends.

Plug-ins outside the style menu (`render.js`, `audio.js`) are classic
scripts. Each is installed by the core with an API that offers:

- the renderer and camera;
- `onScene(fn)` and `onFrame(fn)` hooks;
- the scene data and the clock.

`ctx.rendering` is true during an offline render. A progressive style should
then skip any wait for the view to be still. An option's default `value` may
be a function of the scene (`data`). The ray tracer's accretion disk, for
example, is off by default where the spacetime has no mass.

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
geodesic), `KerrOrbitScene` (`kerr_scene.ndjson`), `TwinParadox` (`twin_paradox.ndjson`),
`BinaryMerger` (`binary_merger.ndjson`; with `neutronStars = true`, two neutron stars), `CosmicExpansion`
(`cosmic_expansion.ndjson`), `CosmicMicrowaveBackground` (`cmb.ndjson`), `GalaxyLensing` (`galaxy_lensing.ndjson`),
`LensInversionLibrary` (`lens_library.ndjson`, eight cases in one gallery),
`LensInversion` (`lens_inversion.ndjson`, plus `lens_inversion_observed.mat`, the synthetic
observation that `LensInversionFromFile` re-analyses into `lens_inversion_file.ndjson`;
`tools/lens_image_to_mat.py` makes such a file from a real picture), `EarthJ2OrbitScene`
(`earth_scene.ndjson`) or `CislunarOrbitScene` (`cislunar_scene.ndjson`). They
are written to `~/Documents/Visualization` when `~/Documents` exists; set the
`sceneDirectory` parameter or the `GR_SCENE_DIR` environment variable to use
another folder, otherwise they go to the simulation working directory.
