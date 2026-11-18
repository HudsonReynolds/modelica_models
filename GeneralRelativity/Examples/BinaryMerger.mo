within GeneralRelativity.Examples;
model BinaryMerger "Two black holes (or neutron stars) spiral together under post-Newtonian gravity with radiation reaction, merge and ring down; recorded for the viewer's ray tracer"
  import GeneralRelativity.CompactBinaries;

  parameter Boolean neutronStars = false "false: binary black hole; true: binary neutron star (finite radius, merger at contact)";
  parameter Real m1Solar = if neutronStars then 1.4 else 36.0 "Mass of body 1 in solar masses (sets Hz, seconds and the neutron-star radius in M)";
  parameter Real m2Solar = if neutronStars then 1.3 else 29.0 "Mass of body 2 in solar masses";
  parameter Real nsRadiusKm = 12.0 "Neutron-star radius" annotation(Dialog(enable = neutronStars));
  parameter Real collapseMassSolar = 2.6 "A neutron-star remnant heavier than this collapses to a black hole" annotation(Dialog(enable = neutronStars));
  parameter Real r0 = if neutronStars then 18.0 else 14.0 "Initial separation in units of M (harmonic coordinates)";
  parameter Real e0(min = 0.0, max = 0.9) = 0.0 "Initial eccentricity: 0 quasi-circular, else start at apocentre";
  parameter Integer pnOrder(min = 0, max = 2) = 2 "Conservative post-Newtonian order of the orbit (radiation reaction is always on)";
  parameter Real rMerge = if neutronStars then R1 + R2 else 3.0
    "Separation at which the bodies merge: contact for neutron stars; about where numerical relativity finds the common horizon for black holes";
  parameter Real sampleInterval = if neutronStars then 1.0 else 0.5 "Time between scene samples, in M";
  parameter String sceneFile = if neutronStars then "binary_neutron_star_merger.ndjson" else "binary_merger.ndjson" "Scene file name";
  parameter String sceneDirectory = "" "Output folder; empty uses GR_SCENE_DIR, else ~/Documents/Visualization, else the working directory";

  final parameter Real M = 1.0 "Total mass: the unit of length and time (G = c = 1)";
  final parameter Real MSolar = m1Solar + m2Solar "Total mass in solar masses";
  final parameter Real m1 = m1Solar / MSolar "Mass of body 1 in M";
  final parameter Real m2 = m2Solar / MSolar "Mass of body 2 in M";
  final parameter Real eta = m1 * m2 "Symmetric mass ratio";
  final parameter Real lengthUnit(unit = "m") = GeneralRelativity.Constants.GM_sun * MSolar / GeneralRelativity.Constants.c^2 "1 M in metres";
  final parameter Real timeUnit(unit = "s") = lengthUnit / GeneralRelativity.Constants.c "1 M in seconds";
  final parameter Real R1 = if neutronStars then nsRadiusKm * 1.0e3 / lengthUnit else 0.0 "Radius of body 1 in M (0: black hole)";
  final parameter Real R2 = if neutronStars then nsRadiusKm * 1.0e3 / lengthUnit else 0.0 "Radius of body 2 in M";
  final parameter Real radiatedFraction = (if neutronStars then 0.5 else 1.0) * CompactBinaries.radiatedEnergyFraction(eta)
    "Mass radiated in gravitational waves (for neutron stars, half the black-hole value: a rough guess)";
  final parameter Real Mf = M * (1.0 - radiatedFraction) "Remnant mass";
  final parameter Real chiF = CompactBinaries.finalSpin(eta) "Remnant spin a/M_f (black-hole fit, also used for a neutron-star remnant)";
  final parameter Real aF = chiF * Mf "Remnant spin length, about +z (the orbital angular momentum)";
  final parameter Boolean remnantIsBlackHole = not neutronStars or MSolar * (1.0 - radiatedFraction) > collapseMassSolar
    "Black holes always merge into one; a neutron-star remnant collapses above collapseMassSolar";
  final parameter Real Rremnant = if remnantIsBlackHole then 0.0 else 1.2 * nsRadiusKm * 1.0e3 / lengthUnit "Radius of a neutron-star remnant in M";
  final parameter Real omegaQNM = CompactBinaries.ringdownOmega(chiF) / Mf "Ringdown angular frequency (black-hole remnant)";
  final parameter Real tauQNM = 2.0 * CompactBinaries.ringdownQuality(chiF) / CompactBinaries.ringdownOmega(chiF) * Mf "Ringdown damping time";
  final parameter Real tAfter = if remnantIsBlackHole then max(12.0 * tauQNM, if neutronStars then 500.0 * M else 0.0) else 100.0 * M
    "Time simulated after the merger: the ringdown; for a neutron-star merger at least 500 M, so the viewer can show the remnant collapsing and the debris settling into a torus";
  final parameter Real tPeters = CompactBinaries.mergerTimePeters(r0, M, eta) * (1.0 - (rMerge / r0)^4)
    "Leading-order (Newtonian + quadrupole) time from r0 to rMerge, for comparison";
  final parameter Real vCircular0 = CompactBinaries.circularOrbitSpeed(r0, M, eta, pnOrder) "Speed of a circular orbit at r0 (exact for the conservative dynamics)";

  Real x[3](start = {r0, 0.0, 0.0}, each fixed = true) "Relative position x1 - x2 (harmonic coordinates, units of M)";
  Real v[3](start = {if e0 > 0.0 then 0.0 else CompactBinaries.inspiralRadialVelocity(r0, M, eta), vCircular0 * sqrt(1.0 - e0), 0.0}, each fixed = true)
    "Relative velocity; counter-clockwise about +z";
  Real tau1(start = 0.0, fixed = true) "Proper time of body 1 (1PN)";
  Real tau2(start = 0.0, fixed = true) "Proper time of body 2 (1PN)";
  Real separation "|x1 - x2|";
  Real omegaOrbit "Orbital angular frequency (x cross v)_z / r^2";
  Real hInspiral[2] "Quadrupole strain {h+, hx} R/M of the inspiral, face-on";
  Real hPlus "Strain h+ times R/M, face-on observer";
  Real hCross "Strain hx times R/M";
  Real fGW(unit = "Hz") "Gravitational-wave frequency: twice the orbital frequency, then the ringdown frequency";
  Real x1[3] "Position of body 1";
  Real x2[3] "Position of body 2";
  discrete Boolean merged(start = false, fixed = true) "The bodies have merged";
  discrete Real tMerge(start = 0.0, fixed = true) "Merger time";
  discrete Real hMergeAmplitude(start = 0.0, fixed = true) "|h| at the merger";
  discrete Real hMergePhase(start = 0.0, fixed = true) "arg(h+ - i hx) at the merger";
  discrete Real omegaPost(start = 0.0, fixed = true) "Angular frequency of the post-merger signal";
  discrete Real tauPost(start = 1.0, fixed = true) "Damping time of the post-merger signal";

  GeneralRelativity.Visualization.SceneRecorder scene(
    spherical = false,
    metricName = "binary",
    p = {m1, m2, R1, R2, 1.0e30, Mf, aF, Rremnant},
    sceneFile = sceneFile,
    sceneDirectory = sceneDirectory,
    title = if neutronStars then "Binary neutron star merger" else "Binary black hole merger",
    referenceRadius = 1.0,
    lengthUnit = "M",
    c = 1.0,
    rInner = 0.5,
    rOuter = r0 + 4.0,
    writeEmbedding = false,
    writeProfile = false,
    writeSlice = false,
    nTracks = 2,
    trackNames = if neutronStars then {"neutron star 1", "neutron star 2"} else {"black hole 1", "black hole 2"},
    trackColors = {"#ffb454", "#5fb3ff"},
    trackDescriptions = {"Heavier body (" + String(m1Solar) + " solar masses); after the merger both tracks sit on the remnant at the origin",
      "Lighter body (" + String(m2Solar) + " solar masses)"},
    nScalars = 4,
    scalarNames = {{"separation", "f_GW Hz", "h+", "hx"}, {"separation", "f_GW Hz", "h+", "hx"}},
    sampleInterval = sampleInterval,
    style = "raytrace",
    trackX = {{time, x1[1], x1[2], x1[3]}, {time, x2[1], x2[2], x2[3]}},
    trackLambda = {tau1, tau2},
    trackScalars = {{separation, fGW, hPlus, hCross}, {separation, fGW, hPlus, hCross}})
    "Scene for the viewer; the metric record {m1, m2, R1, R2, tMerge, Mf, aF, Rremnant} is rewritten at the merger with tMerge";
equation
  separation = sqrt(x * x);
  // The orbit stops at the merger; the remnant sits at the centre of mass (the origin).
  der(x) = if merged then zeros(3) else v;
  der(v) = if merged then zeros(3) else CompactBinaries.relativeAcceleration(x, v, M, eta, pnOrder, true);
  x1 = if merged then zeros(3) else m2 * x;
  x2 = if merged then zeros(3) else -m1 * x;
  omegaOrbit = (x[1] * v[2] - x[2] * v[1]) / separation^2;
  // Clock rates of the bodies at 1PN (own speed and the other body's potential); after the merger the
  // tracks sit inside the remnant and simply follow coordinate time.
  der(tau1) = if merged then 1.0 else 1.0 - 0.5 * (m2 * m2 * (v * v)) - m2 / separation;
  der(tau2) = if merged then 1.0 else 1.0 - 0.5 * (m1 * m1 * (v * v)) - m1 / separation;
  hInspiral = CompactBinaries.quadrupoleStrain(x, v, M, eta);
  // After the merger: a damped sinusoid continuing h = h+ - i hx in amplitude and phase. A black-hole
  // remnant rings at its l = m = 2 quasi-normal mode; a neutron-star remnant is given an illustrative
  // decay at the last inspiral frequency.
  hPlus = if merged then hMergeAmplitude * exp(-(time - tMerge) / tauPost) * cos(hMergePhase - omegaPost * (time - tMerge)) else hInspiral[1];
  hCross = if merged then -hMergeAmplitude * exp(-(time - tMerge) / tauPost) * sin(hMergePhase - omegaPost * (time - tMerge)) else hInspiral[2];
  fGW = if merged then omegaPost / (2.0 * GeneralRelativity.Constants.pi) / timeUnit else omegaOrbit / GeneralRelativity.Constants.pi / timeUnit;
  when separation <= rMerge then
    merged = true;
    tMerge = time;
    hMergeAmplitude = sqrt(hInspiral[1]^2 + hInspiral[2]^2);
    hMergePhase = atan2(-hInspiral[2], hInspiral[1]);
    omegaPost = if remnantIsBlackHole then omegaQNM else 2.0 * omegaOrbit;
    tauPost = if remnantIsBlackHole then tauQNM else 25.0 * M;
  end when;
  assert(separation > 0.5 * rMerge, "The bodies passed through each other without merging");
  assert(merged or time < 50.0 * tPeters + 1000.0, "The binary has not merged");
  when merged and time >= tMerge + tAfter then
    // Within a factor of two of the leading-order estimate; post-Newtonian terms shift the merger time.
    assert(tMerge > 0.5 * tPeters and tMerge < 2.0 * tPeters, "Merger time far from the leading-order (Peters) estimate");
    terminate("Merged and rung down");
  end when;
algorithm
  // The viewer uses the last metric record: rewrite it now that the merger time is known.
  when merged then
    GeneralRelativity.Visualization.writeMetric(GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory), "binary",
      {m1, m2, R1, R2, time, Mf, aF, Rremnant}, "cartesian");
    // Physical scale for the viewer's sound of the chirp: the waveform in real seconds and Hz.
    Modelica.Utilities.Streams.print("{\"type\":\"gwSource\",\"totalMassSolar\":" + GeneralRelativity.Visualization.jsonNumber(MSolar)
      + ",\"secondsPerM\":" + GeneralRelativity.Visualization.jsonNumber(timeUnit)
      + ",\"ringdownHz\":" + GeneralRelativity.Visualization.jsonNumber(if remnantIsBlackHole then omegaQNM / (2.0 * GeneralRelativity.Constants.pi * timeUnit) else 0.0)
      + "}", GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory));
  end when;
  annotation(experiment(StartTime = 0.0, StopTime = 100000.0, Tolerance = 1.0e-10, Interval = 0.5),
    Documentation(info="<html>
<p>A compact binary from inspiral to remnant, in geometric units with the total mass as the unit
(<code>G = c = M = 1</code>; lengths and times in <code>M</code>, which is
<code>GM/c&sup2;</code> = 96 km and 0.32 ms for the default 36 + 29 solar masses, the masses of GW150914).</p>
<ul>
<li><b>Inspiral.</b> The relative orbit <code>x = x<sub>1</sub> &minus; x<sub>2</sub></code> obeys
<a href=\"modelica://GeneralRelativity.CompactBinaries.relativeAcceleration\">CompactBinaries.relativeAcceleration</a>:
Newtonian gravity plus 1PN and 2PN corrections (precession, the faster orbits of strong gravity) in harmonic
coordinates, and the 2.5PN radiation reaction, which makes the orbit shrink as gravitational waves carry
energy away. The start is quasi-circular (the exactly circular speed of the conservative dynamics, <a href=\"modelica://GeneralRelativity.CompactBinaries.circularOrbitSpeed\">circularOrbitSpeed</a>, and the leading-order inward drift), or at
apocentre with eccentricity <code>e0</code>. The bodies sit at <code>x<sub>1</sub> = (m<sub>2</sub>/M)x</code>,
<code>x<sub>2</sub> = &minus;(m<sub>1</sub>/M)x</code>, orbiting counter-clockwise about <code>+z</code>.</li>
<li><b>Waveform.</b> The quadrupole strain <code>h<sub>+</sub>, h<sub>&times;</sub></code> (times
<code>R/M</code>) for an observer on the orbital axis, and the gravitational-wave frequency in Hz: the
chirp.</li>
<li><b>Merger.</b> At separation <code>rMerge</code>: 3 M for black holes (where numerical relativity finds
the common horizon, orbital frequency near 0.18/M), contact <code>R<sub>1</sub> + R<sub>2</sub></code> for
neutron stars. The remnant black hole has the mass and spin of numerical-relativity fits for non-spinning
progenitors (<a href=\"modelica://GeneralRelativity.CompactBinaries.radiatedEnergyFraction\">radiatedEnergyFraction</a>,
<a href=\"modelica://GeneralRelativity.CompactBinaries.finalSpin\">finalSpin</a>): 0.952 M and spin 0.686 for equal
masses. A neutron-star remnant above <code>collapseMassSolar</code> collapses to a black hole; below, it is
drawn as a star of 1.2 times the radius.</li>
<li><b>Ringdown.</b> The strain continues as the remnant's fundamental quasi-normal mode
(<a href=\"modelica://GeneralRelativity.CompactBinaries.ringdownOmega\">ringdownOmega</a>,
<a href=\"modelica://GeneralRelativity.CompactBinaries.ringdownQuality\">ringdownQuality</a>), matched in
amplitude and phase at the merger, for 12 damping times.</li>
</ul>
<p><b>Not modelled:</b> spins of the progenitors, tidal deformation and matter of neutron stars (they are
point masses with a radius; from this record alone the viewer's ray tracer sketches their tidal shapes, the
arms of matter they shed, the collapsing remnant and the debris torus, see the viewer's README), the post-Newtonian breakdown in the last
orbit (the plunge here is only qualitative), the merger itself (the bodies jump to the remnant at the origin)
and recoil.</p>
<p><b>Timing.</b> With the defaults (r0 = 14 M, gravitational-wave frequency 17.5 Hz at the start) the holes
merge after 1708 M (0.55 s), against <code>tPeters</code> = 3030 M for a Newtonian orbit with quadrupole
radiation (<code>pnOrder = 0</code> reproduces that to 0.3%). The 2PN conservative terms flatten the binding
energy toward an innermost stable orbit near <code>r = 6M</code> while the radiation reaction stays at leading
order, so the last orbits shrink faster; a full higher-order inspiral from the same frequency (TaylorT4,
effective-one-body) takes about 3400 M. The equations of motion here are the standard 2.5PN ones of
N-body post-Newtonian codes; use them for the shape of the chirp and the plunge, not for phase-accurate
waveforms. The neutron-star default (1.4 + 1.3 solar masses, 12 km) touches at 6.0 M after 5511 M (73 ms; it starts at a gravitational-wave frequency of 310 Hz)
and its 2.64 solar-mass remnant collapses to a black hole; the scene then continues for 500 M (6.6 ms), long
enough for the viewer to show the remnant collapsing and the debris settling into a torus.</p>
<p>The scene (<code>binary_merger.ndjson</code>, or <code>binary_neutron_star_merger.ndjson</code>) opens in
the viewer's ray tracer, which draws the two holes with the superposed Kerr-Schild metric
<a href=\"modelica://GeneralRelativity.Metrics.binaryKerrSchild\">Metrics.binaryKerrSchild</a> (a snapshot at each
instant: light-travel time across the system is ignored) and the remnant Kerr hole after the merger. Its
metric record is <code>binary</code> with <code>p = {m<sub>1</sub>, m<sub>2</sub>, R<sub>1</sub>, R<sub>2</sub>,
t<sub>merge</sub>, M<sub>f</sub>, a<sub>f</sub>, R<sub>remnant</sub>}</code> (masses in <code>M</code>, radii 0 for black
holes); it is written again at the merger with the merger time.</p>
<p>In the viewer's <b>Embedding / 3D</b> style the scene is the orbital plane as a sheet. It sinks into the
wells of the two masses (by the lapse deficit <code>1 &minus; &alpha;</code> of the superposed metric), and the outgoing
gravitational waves ripple it:</p>
<pre>h(r, &phi;, t) = (M/r) &frac12; Re[(h<sub>+</sub> &minus; i h<sub>&times;</sub>)(t &minus; r/c) e<sup>2i&phi;</sup>]</pre>
<p>This is the leading-order wave in the orbital plane, from the face-on waveform the model records, at the
retarded time. It forms a two-armed spiral moving out at <code>c</code>, chirping as the orbit shrinks, then the
merger and ringdown burst. Its wavelength is about 160 M at the start and 18 M at the merger, so the sheet extends
well beyond the orbit (100 M by default).</p>
</html>"));
end BinaryMerger;
