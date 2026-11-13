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
| `clocks` | Time-contracting grid: a sheet that sinks by the time each static clock has lost, `(1 - alpha)(t - t0)`; ticking clock faces; each body's own proper-time clock | slice |
| `spacetime` | 2+1 spacetime diagram: worldlines against coordinate time; light cones on the "now" plane that narrow and tilt with the metric | tracks; cones need the slice |
| `raytrace` | GPU ray tracer of null geodesics in Kerr-Schild form: shadow, photon ring, lensed sky, thin accretion disk with Doppler and gravitational shifts, ray-traced bodies; orbit or ride along with a body | `metric` record of a known family (`writeMetric`) |

The river, clock and spacetime styles drop one spatial dimension: they draw the
equatorial plane (`theta = pi/2`, or `z = 0` in a Cartesian chart). Weak fields
(Earth, cislunar) have river speeds near `4e-5 c` and clock lags near `1e-9`;
those styles exaggerate them for display and say so in the side panel.

### Style options (ids for `writeView`)

| Style | Option ids (default) |
|---|---|
| `embedding` | `surface` (1), `wire` (1), `fields` (1) |
| `river` | `tracers` (2500), `flow` (1), `flashes` (1), `rings` (5), `map` (1) |
| `clocks` | `warp` (1), `clocks` (1), `stride` (8), `turns` (8), `grid` (1) |
| `spacetime` | `height` (1.6), `cones` (1), `coneSize` (0.05), `rings` (6), `now` (1) |
| `raytrace` | `resolution` (0.6), `steps` (600), `accuracy` (0.04), `disk` (1), `diskOuter` (8, in rs), `shift` (1), `sky` (0 stars, 1 grid, 2 checker), `bodySize` (0.35, in rs), `camera` (0 orbit, k ride with body k), `paths` (0), `exposure` (1) |

Booleans are 0 or 1. For example, a scene that opens in the ray tracer riding
with its first body:

    GeneralRelativity.Visualization.SceneRecorder scene(..., style = "raytrace",
      optionNames = {"camera", "sky"}, optionValues = {1, 2});

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
- optional `update(clock, dt)`, `render(renderer, camera, scene)` (to draw the
  frame itself, as the ray tracer does), `camera()`, `readout()`,
  `onOption(id, value)` (return true when applied in place) and `dispose()`.

`ctx` provides the parsed scene (`data`: header, slice and its fields, tracks,
spheres, metric, embedding rows, profile), the clock, a colour map and colour
bar, a bilinear `sampler(slice, field, mask)` for slice fields, and coordinate
helpers.

## Scene records

`scene` (`beginScene`), `view` (`writeView`), `metric` (`writeMetric`),
`surfaceRow` (`writeEmbeddingSurface`), `profile` (`writeRadialProfile`), `field`
(`writeFieldSlice`), `slice` and `sliceField` (`writeSpacetimeSlice`,
`writeSliceField`), `sphere` (`writeSphere`), `track` (`writeTrack`) and
`sample` (`writeSample`). `Visualization.SceneRecorder` writes all of them from
one component. Samples written without a `track` belong to one default track,
so older scene files still load.

Produce files with `GeneralRelativity.Examples.SchwarzschildOrbitScene`
(`schwarzschild_scene.ndjson`: the Schwarzschild geodesic and the Newton-Cartan
geodesic), `KerrOrbitScene` (`kerr_scene.ndjson`), `EarthJ2OrbitScene`
(`earth_scene.ndjson`) or `CislunarOrbitScene` (`cislunar_scene.ndjson`). They
are written to `~/Documents/Visualization` when `~/Documents` exists; set the
`sceneDirectory` parameter or the `GR_SCENE_DIR` environment variable to use
another folder, otherwise they go to the simulation working directory.
