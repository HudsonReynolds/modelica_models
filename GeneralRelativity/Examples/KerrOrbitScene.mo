within GeneralRelativity.Examples;
model KerrOrbitScene "Equatorial geodesic about a spinning (Kerr) black hole, recorded for every viewer style: frame-dragged river, ZAMO clocks, asymmetric ray-traced shadow"
  extends GeneralRelativity.Geodesics.Geodesic(
    redeclare function metric = GeneralRelativity.Metrics.kerr,
    p = {rs, aSpin},
    c = 1.0,
    x_start = {0.0, r0, GeneralRelativity.Constants.pi / 2.0, 0.0},
    u_start = equatorialTurningPoint(r0, L, rs, aSpin),
    step = {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3},
    massive = true);

  function equatorialTurningPoint "Four-velocity of a massive particle at an equatorial turning point (u^r = 0) with angular momentum L"
    input Real r "Boyer-Lindquist radius";
    input Real L "Specific angular momentum u_phi";
    input Real rs "Schwarzschild radius";
    input Real a "Spin length";
    output Real u[4] "Contravariant four-velocity";
  protected
    Real g[4, 4] = GeneralRelativity.Metrics.kerr({0.0, r, GeneralRelativity.Constants.pi / 2.0, 0.0}, {rs, a});
    Real gInv[4, 4];
    Boolean accepted;
    Real E;
  algorithm
    (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(g);
    // g^tt E^2 - 2 g^tphi E L + g^phiphi L^2 = -1 with u_t = -E, u_phi = L; the root with E > 0.
    E := (gInv[1, 4] * L - sqrt((gInv[1, 4] * L)^2 - gInv[1, 1] * (gInv[4, 4] * L^2 + 1.0))) / gInv[1, 1];
    u := gInv * {-E, 0.0, 0.0, L};
  end equatorialTurningPoint;

  parameter Real M = 1.0 "Black hole mass as a length, G = c = 1";
  final parameter Real rs = 2.0 * M "Schwarzschild radius";
  parameter Real spin = 0.7 "Dimensionless spin a/M, below 1";
  final parameter Real aSpin = spin * M "Spin length a = J/(M c); the name a is the inherited geodesic acceleration";
  final parameter Real rPlus = M + sqrt(M^2 - aSpin^2) "Outer horizon";
  parameter Real r0 = 9.0 * M "Initial radius, a periapsis";
  parameter Real L = 3.6 * M "Specific angular momentum; positive is prograde (with the hole's rotation)";
  final parameter Real E0 = -GeneralRelativity.Geodesics.killingConserved(GeneralRelativity.Metrics.kerr(x_start, p), {1.0, 0.0, 0.0, 0.0}, u_start)
    "Killing energy at the start";
  parameter String sceneFile = "kerr_scene.ndjson" "Scene file name";
  parameter String sceneDirectory = "" "Output folder; empty uses GR_SCENE_DIR, else ~/Documents/Visualization";
  parameter String style = "raytrace" "Style the viewer opens with"
    annotation(choices(choice = "\"embedding\"", choice = "\"river\"", choice = "\"clocks\"", choice = "\"spacetime\"", choice = "\"raytrace\""));
  Real energy "Killing energy -g(d_t, u)";
  Real angularMomentum "Killing angular momentum g(d_phi, u)";
  Real zamoAngularVelocity "Frame-dragging rate -g_tphi/g_phiphi at the particle";

  GeneralRelativity.Visualization.SceneRecorder scene(
    redeclare function metric = GeneralRelativity.Metrics.kerr,
    p = {rs, aSpin},
    step = step,
    spherical = true,
    metricName = "kerr",
    sceneFile = sceneFile,
    sceneDirectory = sceneDirectory,
    title = "Kerr equatorial orbit, a = " + String(spin) + " M",
    referenceRadius = rs,
    lengthUnit = "M",
    c = c,
    rInner = 1.05 * rPlus,
    rOuter = 25.0 * M,
    sliceInner = 0.5 * rPlus,
    sphereNames = {"horizon"},
    sphereRadii = {rPlus},
    sphereColors = {"#000000"},
    nTracks = 1,
    trackNames = {"kerr"},
    trackColors = {"#ff9e64"},
    trackDescriptions = {"Prograde equatorial Kerr geodesic"},
    nScalars = 3,
    scalarNames = {{"E", "L", "omega ZAMO"}},
    sampleInterval = 1.0,
    style = style,
    trackX = {x},
    trackLambda = {time},
    trackScalars = {{energy, angularMomentum, zamoAngularVelocity}});
equation
  energy = -GeneralRelativity.Geodesics.killingConserved(g, {1.0, 0.0, 0.0, 0.0}, u);
  angularMomentum = GeneralRelativity.Geodesics.killingConserved(g, {0.0, 0.0, 0.0, 1.0}, u);
  zamoAngularVelocity = -g[1, 4] / g[4, 4];
  assert(x[2] > rPlus, "The orbit crossed the outer horizon");
  assert(E0 < 1.0, "The orbit is not bound (E >= 1)");
  assert(abs(energy - E0) < 1.0e-6 * E0, "Kerr Killing energy drifted");
  assert(abs(angularMomentum - L) < 1.0e-6 * L, "Kerr Killing angular momentum drifted");
  assert(abs(normResidual) < 1.0e-6, "Four-velocity normalization drifted");
  assert(abs(x[3] - GeneralRelativity.Constants.pi / 2.0) < 1.0e-8, "The orbit left the equatorial plane");
  annotation(experiment(StartTime = 0.0, StopTime = 800.0, Tolerance = 1.0e-10, Interval = 0.5),
    Documentation(info="<html>
<p>A test particle on a bound, prograde, equatorial orbit around a Kerr black hole of spin
<code>a = 0.7M</code>, integrated in Boyer-Lindquist coordinates with the generic
<a href=\"modelica://GeneralRelativity.Geodesics.Geodesic\">Geodesic</a> model. The Killing energy and
angular momentum, the normalization and the equatorial plane are asserted throughout.</p>
<p>The <a href=\"modelica://GeneralRelativity.Visualization.SceneRecorder\">SceneRecorder</a> is the same
component the Schwarzschild example uses, with the metric, its parameters and the ray-tracer family
(<code>kerr</code>) changed. In the viewer:</p>
<ul>
<li><b>raytrace</b> (the default) shows the shadow flattened on the side approaching the camera and
the disk brighter on its approaching side;</li>
<li><b>river</b> shows the flow swirling with the hole: frame dragging, reaching light speed at the
static limit <code>r = 2M</code> on the equator, outside the horizon <code>r<sub>+</sub> &asymp; 1.71M</code>;</li>
<li><b>clocks</b> and <b>spacetime</b> use the clocks and light cones of zero-angular-momentum observers,
whose cones tilt in the direction of rotation.</li>
</ul>
</html>"));
end KerrOrbitScene;
