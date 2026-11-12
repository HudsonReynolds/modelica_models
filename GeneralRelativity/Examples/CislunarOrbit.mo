within GeneralRelativity.Examples;
model CislunarOrbit "High Earth orbit perturbed by the Moon in the weak-field cislunar metric, checked with the Jacobi constant"
  import GeneralRelativity.Constants;
  import GeneralRelativity.Perturbation.Potentials;
  parameter Real GM_earth(unit = "m3/s2") = Constants.GM_earth;
  parameter Real GM_moon(unit = "m3/s2") = Constants.GM_moon;
  parameter Real moonDistance(unit = "m") = Constants.earthMoonDistance;
  parameter Real moonRate(unit = "rad/s") = Constants.omega_moon;
  parameter Real moonPhase0(unit = "rad") = 0.0 "Moon azimuth at t = 0";
  parameter Real c(unit = "m/s") = Constants.c;
  parameter Real r0(unit = "m") = 100.0e6 "Initial circular orbit radius";
  parameter Real phase0(unit = "rad") = Constants.pi / 2.0 "Initial spacecraft azimuth";
  parameter Boolean includeRelativity = true "Add the 1PN Earth term";
  final parameter Real v0(unit = "m/s") = sqrt(GM_earth / r0);
  Real x[3](each unit = "m", start = r0 * {cos(phase0), sin(phase0), 0.0}, each fixed = true) "Position, Earth-centred inertial";
  Real v[3](each unit = "m/s", start = v0 * {-sin(phase0), cos(phase0), 0.0}, each fixed = true) "Velocity";
  Real moon[3](each unit = "m") "Moon position";
  Real r(unit = "m") "Geocentric distance";
  Real moonRange(unit = "m") "Distance to the Moon";
  Real aNewton[3](each unit = "m/s2");
  Real aPostNewtonian[3](each unit = "m/s2");
  Real jacobi(unit = "m2/s2") "Jacobi constant in the frame rotating with the Moon: 2U - |v - Omega z x x|^2";
  Real g[4, 4] "Cislunar metric at the spacecraft";
  Real potentialFromMetric(unit = "m2/s2");
protected
  final parameter Real metricParameters[6] = {GM_earth, GM_moon, moonDistance, moonRate, c, moonPhase0};
  Real potentialParameters[5];
  Real vRotating[3](each unit = "m/s");
  Real effectivePotential(unit = "m2/s2");
  parameter Real jacobi0(unit = "m2/s2", fixed = false) "Jacobi constant at t = 0";
initial equation
  jacobi0 = jacobi;
equation
  moon = Potentials.circularBodyPosition(time, moonDistance, moonRate, moonPhase0);
  potentialParameters = {GM_earth, GM_moon, moon[1], moon[2], moon[3]};
  r = sqrt(x * x);
  moonRange = sqrt((x - moon) * (x - moon));
  aNewton = -Potentials.twoBodyGradient(x, potentialParameters);
  aPostNewtonian = if includeRelativity then GeneralRelativity.Perturbation.postNewtonianAcceleration(x, v, GM_earth, c) else zeros(3);
  der(x) = v;
  der(v) = aNewton + aPostNewtonian;
  vRotating = v - cross({0.0, 0.0, moonRate}, x);
  effectivePotential = 0.5 * moonRate^2 * (x[1]^2 + x[2]^2) - Potentials.twoBody(x, potentialParameters);
  jacobi = 2.0 * effectivePotential - vRotating * vRotating;
  g = GeneralRelativity.Metrics.cislunar({c * time, x[1], x[2], x[3]}, metricParameters);
  potentialFromMetric = GeneralRelativity.Perturbation.newtonianPotential(g, c);
  assert(r > Constants.R_earth and moonRange > 1737.4e3, "The spacecraft collided with a body");
  assert(abs(jacobi - jacobi0) < 1.0e-7 * abs(jacobi0), "Jacobi constant drifted beyond the relativistic correction level");
  // 1 + g_00 is about 1e-10 here, so reading the potential back keeps six digits.
  assert(abs(potentialFromMetric - Potentials.twoBody(x, potentialParameters)) < 1.0e-4 * abs(Potentials.twoBody(x, potentialParameters)),
    "g_00 of Metrics.cislunar does not encode the two-body potential at this time");
  annotation(experiment(StartTime = 0.0, StopTime = 864000.0, Tolerance = 1.0e-10, Interval = 600.0),
    Documentation(info="<html>
<p>The Earth is fixed at the origin and the Moon moves on a circular orbit, so
in the frame rotating with the Moon the Newtonian potential is static and the
Jacobi integral <code>C = 2U<sub>eff</sub> - v<sub>rot</sub>&sup2;</code> is conserved. Its
drift measures the integration error plus the post-Newtonian Earth term, which
is a few parts in <code>10<sup>11</sup></code> at this distance. The metric
<a href=\"modelica://GeneralRelativity.Metrics.cislunar\">Metrics.cislunar</a>
is evaluated along the trajectory and its <code>g<sub>00</sub></code> compared with
the potential. Ten days of flight are simulated; the barycentre motion of the
Earth and the lunar eccentricity are neglected.</p>
<p>This model has no file output, so it builds on every OpenModelica target. <a href=\"modelica://GeneralRelativity.Examples.CislunarOrbitScene\">CislunarOrbitScene</a> adds the Three.js scene stream, which needs <code>Modelica.Utilities</code> to link.</p>
</html>"));
end CislunarOrbit;
