# Three.js viewer for GeneralRelativity scenes

`index.html` renders the newline-delimited JSON scene files written by
`GeneralRelativity.Visualization`. Open it directly in a browser (no server or
build step; Three.js comes from a CDN) and drop a scene file on the page, or use
the file picker. A copy of the viewer can live next to the scene files, e.g. in
`~/Documents/Visualization`:

    cp -R GeneralRelativity/Visualization/threejs/{index.html,render.js,audio.js,render_server.py,render_cli.py,styles} ~/Documents/Visualization/

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
| `raytrace` | GPU ray tracer of null geodesics in Kerr-Schild form: shadow, photon ring, lensed sky, a thin accretion disk glowing as a blackbody (temperature from the hole's mass and accretion rate) with streaks and rings, Doppler and gravitational shifts, ray-traced bodies; orbit or ride along with a body | `metric` record of a known family (`writeMetric`) |
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

When two or more bodies are shown, a legend in the picture's top right corner
lists them: each body's colour and name, plus a short description (the track's
description if short, else its first clause). For example, *schwarzschild*:
Schwarzschild geodesic (general relativity), and *newton*: Newton-Cartan
geodesic (Newtonian gravity). The legend is drawn into the frame, so rendered
videos include it. **legend** in the top bar switches it off. It is left out
in styles that show neither bodies nor tracks; a style can force it with
`frameLegend: true` or `false`.

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
| `raytrace` | `resolution` (1, fraction of device pixels once still), `samples` (24), `steps` (600), `accuracy` (0.04), `disk` (1 around a mass, 0 in flat spacetime), `diskOuter` (8, in rs; in M for a binary, 14 for neutron stars), `diskThick` (0, H/r; above 0 the disk is a volume of gas; 0.22 for neutron stars, where it sets the gas's vertical dispersion), `diskColour` (0 temperature from `diskMass` and `diskRate`, 1 set by `diskTemp`), `diskMass` (log10 solar masses: the scene's `gwSource` mass, else 9.8), `diskRate` (-4.5, log10 of the Eddington fraction), `diskTemp` (5500 K, peak), `diskDetail` (0.8, structure; 0 smooth), `diskMotion` (0; 1 turns the pattern with the clock), `shift` (1), `sky` (0 stars, 1 grid, 2 checker), `bodySize` (0.25, in rs), `bodyLight` (0 the accretion disk, 1 a distant star, 2 the camera), `sunAz` (100°) and `sunEl` (30°, the star's direction), `bodyGlow` (0.06, the bodies' own faint glow), `camera` (0 orbit, k ride with body k), `paths` (0), `tides` (1; neutron stars: tidal shapes, the shed gas as particles that form the torus, the collapsing remnant), `startSep` (a merger scene starts at this separation in M: 12 for neutron stars, else the whole run), `exposure` (1) |
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

### The accretion disk

The disk is a thin disk in the equatorial plane, from the ISCO out to
`diskOuter`, whose gas follows circular geodesics. Each point glows as a
blackbody at the thin-disk temperature (Shakura and Sunyaev, no torque at the
inner edge)

    T(r)^4 = 3 G M Mdot / (8 pi sigma r^3) (1 - sqrt(r_in / r)),

which peaks at `r = 49/36 r_in`. With **disk temperature** set to *from mass
and accretion rate*, the peak follows from the hole's mass (**hole mass**,
log10 M☉) and the accretion rate as a fraction of the Eddington rate
(**accretion**, log10), with `Mdot = L_Edd / (eta c^2)` times that fraction
and the efficiency `eta = 1 - E_isco` of the spin. It scales as
`M^(-1/4) Mdot^(1/4)` per unit Eddington rate: about 5 million K for 10 M☉ at a
tenth of Eddington (X-rays), 30 000 K for 10^8 M☉ (ultraviolet). The default,
M87*'s 6.3 × 10^9 M☉ at 10^-4.5 of Eddington, peaks at 4200 K. Flows that faint
are thought to be hot, thick and optically thin in reality, so take the
default as a thin-disk what-if. Scenes from `Examples.BinaryMerger` start from
the binary's own mass. *Set the peak* takes the temperature directly. The side
panel gives the peak temperature and where its spectrum peaks.

The colour is physical. The specific intensity a camera sees from gas
whose light is shifted by `g` (observed over emitted frequency, from the orbit
and gravity) is a blackbody at `g T`, because `I_nu / nu^3` is invariant along
a ray. The pixel is that blackbody's visible light: Planck's law integrated
against the CIE 1931 colour-matching functions (the fit of Wyman, Sloan and
Shirley, 2013), as linear sRGB. A cool disk is orange, and its approaching side
is far brighter and whiter than the receding side. A disk hotter than about
30 000 K shows only the blue-white Rayleigh-Jeans tail in visible light, where
the brightness follows `g T` and the colour hardly changes. The brightness unit
is the disk's mean visible luminance out to 4 `r_in`, so hot and cool disks
show at a similar overall brightness, like a camera's automatic exposure, while
the variation across the disk stays physical. **exposure** scales it.

**disk structure** adds texture that is not computed physics, so the disk reads
as a disk:
- streaks along the orbits, about ten per turn and twenty per e-fold of
  radius;
- the trailing-spiral shear of differential rotation (each radius turns at its
  orbital angular velocity; the pattern shows four radians of turn at the inner
  edge);
- rings and gaps in radius;
- warm and cool patches of ±12% in temperature.

It sets the fraction of the surface that glows and blocks light, so the sky
and the lensed far side show through the gaps. The pattern is fixed. With
**disk rotates**, or the same switch in the Render dialog, it turns with the
clock at each radius's own speed. To keep the winding bounded, two copies of
the pattern, each living one inner orbit, cross-fade with their contrast kept
steady. Every clock step then re-traces the whole image, so this is meant for
rendered videos.

**disk thickness** (H/r) above 0 turns the disk into a volume of gas: a
Gaussian in height with scale height `H = (H/r) r`, with an opacity that gives a
vertical optical depth of `3 r_in / r`. It is opaque only near the hole, so
elsewhere its brightness follows the density of its structure: a glowing cloud.
It glows as a blackbody at the same temperature profile, Doppler shifted by
circular orbits about the axis, with the midplane up to 1.25 times hotter than
the cooler atmosphere. Its source function is 0.22 of the blackbody's
brightness, against an opaque thin disk's 0.4, because a thick torus shows
much more glowing surface. Each integration step adds emission and absorption
at one point placed at random along the step (re-drawn every frame); the
refinement samples average the placements. The structure is sheared by the
differential rotation like the disk's: a two-armed spiral, which survives the
long sight lines through the gas, plus streaks and rings, which mostly average
out. A thick disk takes up to twice as long to trace as a thin one.

### Lighting the bodies

The bodies are lit along straight light paths (not lensed, which would cost a
second ray trace per pixel). They show their track colour as albedo times the
light, plus a faint glow of their own (**night-side glow**, not physical), so
the dark side stays visible. **bodies lit by** picks the light:
- **The accretion disk** (the default). The disk acts as 6 rings of 16 glowing
  points, each with the disk's blackbody colour and brightness at its radius
  and the mean coverage of its structure. Light is blocked where its straight
  path passes within 1.5 horizon radii of the hole, about the photon sphere.
  The points shine in all directions, because a thin disk seen edge-on would
  leave a body in its own plane dark. The side of a body facing the hot inner
  disk is brightest. The light is tone-mapped (`1 - exp(-E)`), so a body far
  from the disk stays visible while keeping the light's colour and direction.
  Without a disk, a distant star lights the bodies instead.
- **A distant star** in the direction **star azimuth**, **star elevation**
  (degrees, in the scene's frame, z along the spin). The hole can eclipse it,
  and where the path crosses the disk, the disk's structure shades it.
- **The camera**, the earlier flat shading.

### The sky

The ray tracer's sky has these parts:
- **Stars:** two layers, a sparse one and a dense one that is denser towards a
  Milky Way band. Brightness follows `N(>F) ~ F^-1.5` (many faint, few bright),
  and colours are blackbodies of 3000 to 30 000 K.
- **Grid** (sky option 1): lines of latitude and longitude every 15°.
- **Checker** (sky option 2).

Stars are point sources, as in the renderer of James et al. (Class. Quantum
Grav. 32, 065001, 2015). Each pixel's footprint on the sky comes from
screen-space derivatives of the rays' final directions:
- every star is drawn as a dot about a pixel across wherever lensing puts it,
  rather than being stretched into an arc near the Einstein ring;
- its brightness is multiplied by the magnification, the solid angle of the
  pixel at the camera over that on the sky, capped at 20;
- grid lines are about 1.5 pixels wide at any resolution, since a line's
  distance is divided by its change per pixel across it.

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

With **tides and debris** on (the default), a neutron-star binary (`R_a > 0`)
is drawn with matter, as far as a ray tracer without hydrodynamics can. Every
shape, level and density below changes continuously with the clock, so the
phases follow each other without jumps:
- **Tidal shapes.** The stars are one surface, below an equipotential of the
  Roche potential: the two point masses in the frame that turns with the
  recorded orbit about the centre of mass. Each star's level is the potential
  at its radius on the side away from its companion, so the side facing the
  companion bulges towards the inner Lagrange point L1 (a teardrop). Before
  either level reaches L1's, the two blend into their mean. The stars then
  reach L1 together, and the neck between them opens and widens into a peanut
  (from about 7.9 M) until the model declares the merger at 6 M. Matter at the
  neck is drawn hotter and whiter.
- **Shedding.** Once the surfaces come within 12% of L1's potential (at a
  separation of 9.7 M, 1.5 orbits and 250 M before the merger), the stars'
  outer tips shed gas, 60% of it from the lighter star. The gas is 30 000 test
  particles (see *The gas* below), launched at a rate that rises until the
  merger and stops as the remnant closes up. Each leaves a tip, spread over its
  face, moving with the turning pair plus an outward kick of 0.15 to 0.35 c.
  The tips turn faster than the gas, so the streams trail as spiral arms.
  Freshly shed gas has its star's colour and turns to the gas's blackbody
  colour over 0.8 orbits.
- **The remnant.** After the merger the pair keeps turning at the orbit's rate
  and closes up over 0.6 orbits: its separation continues smoothly from the
  plunge's value and rate to zero, and its shape blends into a ball. Meanwhile
  the metric changes from the two-centre one, with the masses riding the lobes,
  to the remnant's Kerr metric (`Mf`, `aF`): the two masses' terms fade out as
  the Kerr-Schild term of the spinning hole fades in (`g = eta + (1 - w)(f1 k1
  k1 + f2 k2 k2) + w f k k`, again a Kerr-Schild superposition). The lensed sky
  and grid therefore slow from turning with the pair to the spinning hole's
  fixed pattern in one motion, without first stopping as a non-spinning hole.
  The ball then shrinks onto the horizon over 0.9 orbits.
- **The collapse.** The ball is dimmed by its gravitational redshift: the cube
  of `sqrt(-g_tt)` relative to the camera's, and to 0.7, about a neutron star's
  own. It fades to black as it reaches the horizon, so the shadow appears
  without a jump. Inside the photon sphere it is seen at the shadow's size.
- **The gas.** The particles move under gravity alone: the pair as two
  softened point masses on the recorded orbit, then closing up, then the
  remnant's mass at the centre. After the merger they dissipate as shocks
  would. Each particle's radial velocity is damped within about half an
  orbit while its angular momentum is kept, so the gas settles where its
  angular momentum puts it (mostly 5 to 14 M). Its vertical motion relaxes to
  a dispersion that holds it at **disk thickness** (H/r, 0.22 by default),
  and a little in-plane noise spreads it. Gas that falls back into a star or
  the remnant, or inside the hole's innermost stable orbit, is removed: about
  four fifths of it, leaving some 6000 particles in the torus. So the arms
  leaving the stars are the gas that wraps around the remnant and becomes
  the torus. Nothing appears from nowhere, and the torus rotates because its
  gas orbits.
- **Drawing it.** Each frame the particles are spread over a 128 x 128 x 32
  grid around the hole (a Gaussian about 0.5 M wide). The grid holds the
  density, plus the density of young gas from each star, for its colour. The
  ray tracer samples it as a 3D texture, with the disk's finer streaks on top
  at half strength. It is drawn out to **disk outer** (14 M by default),
  fading beyond. The gas emits as a blackbody at the thin-disk temperature
  for its radius, with the Doppler shift of circular orbits. Its temperature
  comes from the remnant's mass by default (blue-white at any realistic
  accretion rate); *set the peak* gives a cooler colour.
- **Cost.** The simulation uses leapfrog steps of 0.3 M and is deterministic,
  with a checkpoint every 19 M so the clock can go back. Jumping straight to
  the end of the run takes about 1.6 s once; after that, playing or rendering
  forward costs about 15 ms per frame at 4 M per frame, and stepping back
  about as much. The side panel counts the particles shed, in flight and
  removed.

None of it is a fluid simulation: the shapes are Newtonian equipotentials
and the gas is test particles with a dissipation rule, a sketch of what
merger simulations show.
`Examples.BinaryMerger` records 500 M after a neutron-star merger so the torus
has time to form. The disk is suppressed before the merger, so a black-hole
binary with the disk on gets its disk at the merger. The side panel names the
phase: inspiral, shedding or contact, closing up, collapsing (with its
progress) and the settling debris.

A merger scene starts at **start at separation**: 12 M for neutron stars,
3.5 orbits and 750 M before the merger, instead of the recorded 18 M (5511 M
before it); the whole run at or above the first recorded separation. The time
slider and the Render dialog use that range.

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
- **Clock:** steps through a time range (by default the viewer's, which for a
  merger begins at **start at separation**) at a frame rate (24 fps by
  default) and a video length that you choose.
- **Slow motion** (merger scenes): frames are this many times denser in
  simulation time around the merger, easing in and out over half that span on
  each side. For neutron stars the span runs from the first shedding through
  the collapse and the torus's filling (three orbits after the merger), and
  the default is 4×. For black holes it runs from the last 1.5 orbits to 100 M
  of ringdown, and the default is 1× (off). The dialog shows the resulting
  simulation time per frame. 1 keeps the clock steady. Live playback is not
  slowed.
- **Waiting:** each frame waits until the active style reports it final. For
  the ray tracer, that means every refinement sample (**samples per frame**,
  8 by default) plus the bodies' layer. The GPU does this work and is fully
  busy while rendering (98 to 100% on an Apple M4), so time per frame is set by
  the GPU and by the samples: about 2.5 s per 1080p frame of the neutron-star
  merger at 6 samples, 4 s at 12 and 7.5 s at 24. Against 48 samples, 6 differ
  by 0.19 of 255 grey levels on average and 12 by 0.13, far below what video
  compression discards, so more samples only pay for stills.
- **Resolution:** set in the dialog, independent of the window. The view shows
  the frame letterboxed while it renders.
- **Legend:** the in-picture legend of the bodies is part of every frame
  while **legend** is ticked in the top bar.
- **Accretion disk:** the ray tracer's disk can turn with the clock for the
  render (**the accretion disk rotates**) and stays fixed in the viewer.

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

### From a terminal

`render_cli.py` renders without opening the viewer, with the same renderer
and the same frames as the dialog (pixel for pixel). It serves the viewer
folder, opens it in a headless browser (Firefox, else Chrome or Chromium),
loads the scene and settings, and writes the frames and the video through the
helper. The render keeps going with the browser closed and does not pause
when other windows are in front.

    cd ~/Documents/Visualization
    python3 render_cli.py binary_neutron_star_merger.ndjson -o renders/merger.mp4
    python3 render_cli.py kerr_scene.ndjson --style raytrace --seconds 20 --size 3840x2160 -o kerr.mp4

The dialog's **command line** row (with **Copy**) gives the command for what
is set up in the viewer: the style, the options that differ from their
defaults, hidden bodies, the legend, the camera, the time range, frames and
output. Settings (`python3 render_cli.py --help` lists them):

- `--style`, `--option KEY=VALUE` (repeatable), `--hide NAME`, `--no-legend`;
- `--camera x,y,z` and `--target x,y,z` (viewer coordinates; default: the
  style's camera);
- `--from`, `--to` (clock; default: the viewer's range), `--fps` (24),
  `--seconds` (10), `--size` (1920x1080), `--samples` (8), `--slow` (the
  style's), `--disk-rotates`;
- `-o`: `.mp4` (H.264), `.mpeg` (MPEG-2), or a folder ending in `/` for PNG
  frames; `--keep-frames`.

Frames are independent and reproducible, including a merger's gas, so a long
render can be split across machines (a faster GPU renders faster; the
`WebGL:` line names the renderer the browser uses, and a software one such as
SwiftShader or llvmpipe is very slow). Each part writes PNG frames numbered
within the whole sequence, and `--encode` joins them:

    python3 render_cli.py scene.ndjson <same settings> --frames 0:120 -o frames/      # machine 1
    python3 render_cli.py scene.ndjson <same settings> --frames 120:240 -o frames/    # machine 2
    python3 render_cli.py --encode frames/ --fps 24 -o merger.mp4

Pages can render the same way through `GRRender.render(settings)` (see
`render.js`), which `render_cli.py` calls.

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
example, is off by default where the spacetime has no mass. An optional
`show(data, values)` hides the option when it returns false. `values` holds the
style's current option values, and the panel is rebuilt when a checkbox or a
choice changes; the disk's own options, for example, are listed only while the
disk is on.

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
