within GeneralRelativity.Examples;
model EarthJ2Orbit "Low Earth orbit in the weak-field Earth metric: J2, post-Newtonian and Lense-Thirring accelerations with the tidal field and a scene stream"
  import GeneralRelativity.Constants;
  import GeneralRelativity.Perturbation.Potentials;
  parameter Real GM(unit = "m3/s2") = Constants.GM_earth "Gravitational parameter";
  parameter Real R(unit = "m") = Constants.R_earth "Equatorial radius";
  parameter Real J2 = Constants.J2_earth "Second zonal harmonic";
  parameter Real J(unit = "kg.m2/s") = Constants.J_earth "Spin angular momentum about +z";
  parameter Real c(unit = "m/s") = Constants.c "Speed of light";
  parameter Real semiMajorAxis(unit = "m") = R + 700.0e3 "Semi-major axis";
  parameter Real eccentricity = 0.01;
  parameter Real inclination(unit = "rad") = 51.6 * Constants.pi / 180.0 "Inclination to the equator";
  parameter Boolean includeRelativity = true "Add the 1PN and Lense-Thirring accelerations";
  parameter String sceneFile = "earth_scene.ndjson" "Scene stream for Visualization/threejs/index.html";
  parameter Real sampleInterval(unit = "s") = 60.0;
  final parameter Real p[3] = {GM, R, J2} "Potential parameters";
  final parameter Real perigee(unit = "m") = semiMajorAxis * (1.0 - eccentricity);
  final parameter Real perigeeSpeed(unit = "m/s") = sqrt(GM * (1.0 + eccentricity) / perigee);
  final parameter Real meanMotion(unit = "rad/s") = sqrt(GM / semiMajorAxis^3);
  final parameter Real nodalRegressionRate(unit = "rad/s") =
    -1.5 * meanMotion * J2 * (R / (semiMajorAxis * (1.0 - eccentricity^2)))^2 * cos(inclination)
    "Secular J2 node rate, first order";
  final parameter Real relativisticPerigeeAdvance(unit = "rad") =
    6.0 * Constants.pi * GM / (c^2 * semiMajorAxis * (1.0 - eccentricity^2)) "GR perigee advance per orbit";
  Real x[3](each unit = "m", start = {perigee, 0.0, 0.0}, each fixed = true) "Position, Earth-centred inertial";
  Real v[3](each unit = "m/s", start = {0.0, perigeeSpeed * cos(inclination), perigeeSpeed * sin(inclination)}, each fixed = true)
    "Velocity";
  Real r(unit = "m") "Geocentric distance";
  Real aNewton[3](each unit = "m/s2") "Point mass plus J2";
  Real aPostNewtonian[3](each unit = "m/s2");
  Real aLenseThirring[3](each unit = "m/s2");
  Real relativeCorrection "Size of the relativistic acceleration relative to the Newtonian one";
  Real newtonianEnergy(unit = "m2/s2") "v^2/2 + Phi; conserved by the static Newtonian part";
  Real h[3](each unit = "m2/s") "Specific angular momentum";
  Real raan(unit = "rad") "Right ascension of the ascending node";
  Real g[4, 4] "Weak-field Earth metric at the satellite";
  Real potentialFromMetric(unit = "m2/s2") "Phi read back from g_00";
  Real tidal[3, 3](each unit = "1/s2") "Newtonian tidal tensor at the satellite";
  Real tidalLaplacian(unit = "1/s2") "Trace of the tidal tensor, zero in vacuum";
  Real einsteinLinear[4, 4] "Linearized Einstein tensor at the satellite, zero in vacuum";
protected
  final parameter Real metricParameters[5] = {GM, R, J2, c, J};
  final parameter Real newtonianEnergy0(unit = "m2/s2") = 0.5 * perigeeSpeed^2 + Potentials.oblateJ2({perigee, 0.0, 0.0}, p);
  final parameter Real sliceCoordinate[40] = {-3.0 * R + 6.0 * R * (i - 1) / 39.0 for i in 1:40} "Even count so no sample sits on the singular origin";
  Real riemannLinear[4, 4, 4, 4];
initial algorithm
  GeneralRelativity.Visualization.beginScene(sceneFile, "Earth J2 orbit", "cartesian", R, "m", c);
  GeneralRelativity.Visualization.writeSphere(sceneFile, "Earth", R, zeros(3), "#3b7dd8");
  GeneralRelativity.Visualization.writeFieldSlice(sceneFile, "tidal magnitude |Phi_,ij|", "xy", 0.0,
    sliceCoordinate, sliceCoordinate,
    GeneralRelativity.Visualization.tidalMagnitudeSlice(Potentials.oblateJ2Gradient, p, 1.0, sliceCoordinate, sliceCoordinate, 0.0, R));
equation
  r = sqrt(x * x);
  aNewton = -Potentials.oblateJ2Gradient(x, p);
  aPostNewtonian = if includeRelativity then GeneralRelativity.Perturbation.postNewtonianAcceleration(x, v, GM, c) else zeros(3);
  aLenseThirring = if includeRelativity then GeneralRelativity.Perturbation.lenseThirringAcceleration(x, v, {0.0, 0.0, J}, c) else zeros(3);
  der(x) = v;
  der(v) = aNewton + aPostNewtonian + aLenseThirring;
  relativeCorrection = sqrt((aPostNewtonian + aLenseThirring) * (aPostNewtonian + aLenseThirring)) / sqrt(aNewton * aNewton);
  newtonianEnergy = 0.5 * v * v + Potentials.oblateJ2(x, p);
  h = cross(x, v);
  raan = atan2(h[1], -h[2]);
  g = GeneralRelativity.Metrics.earth({c * time, x[1], x[2], x[3]}, metricParameters);
  potentialFromMetric = GeneralRelativity.Perturbation.newtonianPotential(g, c);
  tidal = Potentials.hessianFromGradient(Potentials.oblateJ2Gradient, x, p, 1.0);
  (riemannLinear, einsteinLinear, tidalLaplacian) = GeneralRelativity.Perturbation.weakFieldCurvature(tidal, c);
  assert(r > R, "The satellite hit the Earth");
  assert(abs(newtonianEnergy - newtonianEnergy0) < 1.0e-7 * abs(newtonianEnergy0),
    "Newtonian energy drifted beyond the size of the relativistic corrections");
  assert(not includeRelativity or (relativeCorrection > 1.0e-10 and relativeCorrection < 1.0e-8),
    "Relativistic corrections are outside the expected 1e-9 range for a low Earth orbit");
  // 1 + g_00 is about 1.4e-9 here, so reading the potential back keeps seven digits.
  assert(abs(potentialFromMetric - Potentials.oblateJ2(x, p)) < 1.0e-5 * abs(Potentials.oblateJ2(x, p)),
    "g_00 of Metrics.earth does not encode the J2 potential");
  assert(abs(tidalLaplacian) < 1.0e-6 * GM / r^3, "The J2 potential is not harmonic outside the Earth");
algorithm
  when sample(0.0, sampleInterval) then
    GeneralRelativity.Visualization.writeSample(sceneFile, time, {c * time, x[1], x[2], x[3]},
      {relativeCorrection, raan, tidalLaplacian});
  end when;
  annotation(experiment(StartTime = 0.0, StopTime = 12000.0, Tolerance = 1.0e-10, Interval = 10.0),
    Documentation(info="<html>
<p>The equations of motion are the Newtonian gradient of the J2 potential plus
the two first-order relativistic terms that follow from the weak-field metric
<a href=\"modelica://GeneralRelativity.Metrics.earth\">Metrics.earth</a>: the
Schwarzschild (post-Newtonian) term and the Lense-Thirring term from
<code>g<sub>0i</sub></code>. Their combined size is about <code>10<sup>-9</sup></code>
of the Newtonian acceleration. The node rate from J2 is about
<code>-5 degrees per day</code> for this orbit; the relativistic perigee advance is
<code>relativisticPerigeeAdvance</code> per revolution, some 0.01 arcseconds.</p>
<p>At every step the metric is evaluated at the satellite and its
<code>g<sub>00</sub></code> read back as a potential; the Hessian of the potential
gives the tidal tensor and the linearized Einstein tensor, whose vanishing
confirms the J2 field is a vacuum solution of the linearized equations.</p>
<p>The scene stream contains the Earth, the orbit sampled every minute and the
tidal magnitude on the equatorial plane out to three Earth radii.</p>
</html>"));
end EarthJ2Orbit;
