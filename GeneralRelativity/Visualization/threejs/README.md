# Three.js viewer for GeneralRelativity scenes

`index.html` is self-contained (the viewer script is inline) and renders the newline-delimited JSON scene files written by
`GeneralRelativity.Visualization` (`beginScene`, `writeEmbeddingSurface`,
`writeRadialProfile`, `writeFieldSlice`, `writeSphere`, `writeSample`).

Open `index.html` directly in a browser and drop a scene file on the page, or
use the file picker. Three.js is loaded from a CDN by the import map in
`index.html`; no build step or local server is required.

What it shows:

- the embedding surface of the equatorial plane (Flamm's paraboloid for
  Schwarzschild), coloured by the logarithm of the Kretschmann scalar;
- reference spheres such as the horizon, the photon sphere or a planet;
- the integrated trajectory as a line, with a particle animated against the
  affine parameter;
- a scalar field slice (for example the tidal magnitude near Earth) as a
  coloured plane;
- a side panel plotting `g_tt`, `g_rr`, the Ricci scalar and the Kretschmann
  scalar against radius.

Spherical samples `{ct, r, theta, phi}` are mapped to Cartesian coordinates and
lifted onto the embedding surface when one is present; Cartesian samples are
drawn directly.

Produce a file with `GeneralRelativity.Examples.SchwarzschildOrbitScene`
(`schwarzschild_scene.ndjson` in the simulation working directory),
`EarthJ2OrbitScene` (`earth_scene.ndjson`) or `CislunarOrbitScene`
(`cislunar_scene.ndjson`).
