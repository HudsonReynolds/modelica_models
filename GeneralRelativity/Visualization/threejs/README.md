# Three.js viewer for GeneralRelativity scenes

`index.html` is self-contained (the viewer script is inline) and renders the newline-delimited JSON scene files written by
`GeneralRelativity.Visualization` (`beginScene`, `writeEmbeddingSurface`,
`writeRadialProfile`, `writeFieldSlice`, `writeSphere`, `writeTrack`, `writeSample`).

Open `index.html` directly in a browser and drop a scene file on the page, or
use the file picker. Three.js is loaded from a CDN by the import map in
`index.html`; no build step or local server is required.

What it shows:

- the embedding surface of the equatorial plane (Flamm's paraboloid for
  Schwarzschild), coloured by the logarithm of the Kretschmann scalar;
- reference spheres such as the horizon, the photon sphere or a planet;
- every trajectory track as an animated body with a growing trail, all moving
  together against coordinate time `t = x[0]/c` (so a proper-time geodesic and a
  Newtonian body sampled at the same `ct` stay in step), the full path drawn
  faintly, and a dashed link between the first two bodies. Playback starts as
  soon as a file is loaded; the `t` slider scrubs, `speed` sets the rate (the
  whole run takes about a minute at the default), and `loop` repeats it;
- a legend with each body's colour, description and live readout (affine
  parameter, position and the named scalars declared by `writeTrack`);
- a scalar field slice (for example the tidal magnitude near Earth) as a
  coloured plane;
- a side panel plotting `g_tt`, `g_rr`, the Ricci scalar and the Kretschmann
  scalar against radius.

Spherical samples `{ct, r, theta, phi}` are mapped to Cartesian coordinates and
lifted onto the embedding surface when one is present; Cartesian samples are
drawn directly.

Produce a file with `GeneralRelativity.Examples.SchwarzschildOrbitScene`
(`schwarzschild_scene.ndjson`: the Schwarzschild geodesic and the Newton-Cartan
geodesic from the same initial data),
`EarthJ2OrbitScene` (`earth_scene.ndjson`) or `CislunarOrbitScene`
(`cislunar_scene.ndjson`).

Scene files are written to `~/Documents/Visualization` when `~/Documents`
exists. Set the `sceneDirectory` parameter or the `GR_SCENE_DIR` environment
variable to use another folder; otherwise they go to the simulation working
directory.

Samples written without a `track` belong to one default track, so older scene
files still load.
