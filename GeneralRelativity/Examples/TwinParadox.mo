within GeneralRelativity.Examples;
model TwinParadox "Twin paradox in Minkowski space: one twin stays home, the other flies out and back under constant proper acceleration; ages integrated, checked against hyperbolic motion and recorded for the interactive diagram"
  function atanh_ "Inverse hyperbolic tangent"
    input Real x;
    output Real y;
  algorithm
    y := 0.5 * log((1.0 + x) / (1.0 - x));
  end atanh_;

  function acosh_ "Inverse hyperbolic cosine, x >= 1"
    input Real x;
    output Real y;
  algorithm
    y := log(x + sqrt(x^2 - 1.0));
  end acosh_;

  parameter Real beta(min = 0.0, max = 0.9999) = 0.9 "Cruise speed v/c";
  parameter Real distance(unit = "ly") = 4.0 "Turnaround distance from home";
  parameter Real accelerationG = 1.0 "Proper acceleration felt by the traveller, in g";
  constant Real gInLyPerYr2 = 9.80665 * (365.25 * 86400.0)^2 / 9.4607304725808e15 "1 g in light-years per year squared (about 1.032)";
  final parameter Real a = accelerationG * gInLyPerYr2 "Proper acceleration, ly/yr^2 (c = 1)";
  final parameter Real etaMax = min(atanh_(beta), acosh_(1.0 + a * distance / 2.0))
    "Peak rapidity: the cruise rapidity, or less when the turnaround comes before cruise speed is reached";
  final parameter Real tauAccelerate = etaMax / a "Traveller's proper time for one acceleration leg";
  final parameter Real xAccelerate = (cosh(etaMax) - 1.0) / a "Distance covered in one acceleration leg";
  final parameter Real tauCruise = max(0.0, (distance - 2.0 * xAccelerate) / tanh(etaMax)) / cosh(etaMax)
    "Traveller's proper time for one cruise leg";
  final parameter Real tauReunion = 4.0 * tauAccelerate + 2.0 * tauCruise "Traveller's age at the reunion (closed form)";
  final parameter Real tReunion = 4.0 * sinh(etaMax) / a + 2.0 * tauCruise * cosh(etaMax) "Home twin's age at the reunion (closed form)";

  Real x[4](start = {0.0, 0.0, 0.0, 0.0}, each fixed = true) "Traveller's event {ct, x, y, z} in light-years; Modelica time is the traveller's proper time in years";
  Real u[4](start = {1.0, 0.0, 0.0, 0.0}, each fixed = true) "Traveller's four-velocity dx/d tau";
  Real alpha "Proper acceleration along +x, the reading of the traveller's accelerometer";
  Real normResidual "g(u,u) + 1, zero for a correctly normalized four-velocity";
  Real homeAge "Home twin's age: their proper time is the inertial coordinate time x[1]/c";
  Real travellerAge "Traveller's age: Modelica time";
  Real ageDifference "homeAge - travellerAge";
  Real speed "Traveller's coordinate speed dx/dt in units of c";

  GeneralRelativity.Visualization.SceneRecorder scene(
    redeclare function metric = GeneralRelativity.Metrics.minkowski,
    p = {0.0},
    spherical = false,
    metricName = "minkowski",
    sceneFile = "twin_paradox.ndjson",
    title = "Twin paradox",
    referenceRadius = 1.0,
    lengthUnit = "ly",
    c = 1.0,
    rInner = 0.05,
    rOuter = distance + 1.0,
    writeEmbedding = false,
    writeProfile = false,
    writeSlice = false,
    nTracks = 2,
    trackNames = {"home", "traveller"},
    trackColors = {"#5fb3ff", "#ffb454"},
    trackDescriptions = {"Stays home: an inertial worldline, the longest proper time between the two meetings",
      "Flies out and back with constant proper acceleration"},
    nScalars = 1,
    scalarNames = {{"age"}, {"age"}},
    sampleInterval = 0.02,
    style = "twins",
    optionNames = {"beta", "distance", "acceleration"},
    optionValues = {beta, distance, accelerationG},
    trackX = {{x[1], 0.0, 0.0, 0.0}, x},
    trackLambda = {homeAge, travellerAge},
    trackScalars = {{homeAge}, {travellerAge}}) "Scene for the viewer's twin-paradox diagram";
equation
  // Out at +a, cruise, -a for twice as long (stop, then head home), cruise, +a to arrive at rest.
  alpha = if time < tauAccelerate then a
    else if time < tauAccelerate + tauCruise then 0.0
    else if time < 3.0 * tauAccelerate + tauCruise then -a
    else if time < 3.0 * tauAccelerate + 2.0 * tauCruise then 0.0
    else if time < tauReunion then a
    else 0.0;
  der(x) = u;
  // Four-acceleration alpha (u^x, u^t, 0, 0): orthogonal to u, with magnitude alpha.
  der(u) = alpha * {u[2], u[1], 0.0, 0.0};
  normResidual = -u[1]^2 + u[2]^2 + u[3]^2 + u[4]^2 + 1.0;
  homeAge = x[1];
  travellerAge = time;
  ageDifference = homeAge - travellerAge;
  speed = u[2] / u[1];
  assert(abs(normResidual) < 1.0e-6, "Traveller's four-velocity lost its normalization");
  assert(x[2] <= distance + 1.0e-6, "The traveller went beyond the turnaround distance");
  when time >= tauReunion then
    assert(abs(x[2]) < 1.0e-6, "The traveller did not return home");
    assert(abs(u[1] - 1.0) < 1.0e-6 and abs(u[2]) < 1.0e-6, "The traveller did not arrive at rest");
    assert(abs(homeAge - tReunion) < 1.0e-6 * tReunion, "Home twin's age at the reunion disagrees with hyperbolic motion");
    assert(ageDifference > 0.0, "The travelling twin is not the younger one");
    terminate("The twins have met again");
  end when;
  annotation(experiment(StartTime = 0.0, StopTime = 200.0, Tolerance = 1.0e-10, Interval = 0.01),
    Documentation(info="<html>
<p>Flat spacetime, inertial coordinates <code>{ct, x, y, z}</code> in years and light-years
(<code>c = 1</code>). One twin stays at <code>x = 0</code>; the other flies to <code>x = D</code> and back,
with an accelerometer that reads <code>+a</code>, <code>0</code>, <code>&minus;a</code>, <code>0</code>,
<code>+a</code> in five legs. Modelica <code>time</code> is the traveller's proper time; the model integrates
<code>dx/d&tau; = u</code> and <code>du/d&tau; = &alpha;(u<sup>x</sup>, u<sup>t</sup>, 0, 0)</code>, a
four-acceleration orthogonal to <code>u</code> (hyperbolic motion during each acceleration leg).</p>
<p>The home twin's proper time is the coordinate time. At the reunion the model asserts the closed forms
<code>t = 4 sinh(&eta;)/a + 2&tau;<sub>c</sub> cosh &eta;</code> and <code>&tau; = 4&eta;/a + 2&tau;<sub>c</sub></code>
(rapidity <code>&eta; = atanh(v/c)</code>, cruise proper time <code>&tau;<sub>c</sub></code>), that the traveller is
home and at rest, and that the traveller is younger. With the defaults (0.8 c, 4 ly, 1 g) the home twin
is 11.94 years older at the reunion and the traveller 8.32 years.</p>
<p>The traveller's worldline is the shorter proper time because it is not a geodesic: of all timelike paths
between the two meetings, the inertial (geodesic) one has the longest proper time. Nothing is symmetric:
only the traveller feels acceleration.</p>
<p>The scene opens in the viewer's <b>twin paradox</b> style: a Minkowski diagram whose sliders (speed,
distance, acceleration) recompute the worldlines live with the same closed forms, with proper-time ticks,
yearly light signals and the traveller's line of simultaneity, which swings across the home worldline
during the turnaround. The simulated worldlines are overlaid for comparison.</p>
</html>"));
end TwinParadox;
