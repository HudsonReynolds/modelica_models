within GeneralRelativity.Examples;
model PolytropeStar "Tolman-Oppenheimer-Volkoff integration of a polytropic star, checked against the n = 1 Lane-Emden solution; Modelica time is the radius"
  import pi = GeneralRelativity.Constants.pi;
  parameter Real K = 1.0 "Polytropic constant in geometric units (G = c = 1)";
  parameter Real Gamma = 2.0 "Adiabatic index; the Lane-Emden check needs Gamma = 2 (n = 1)";
  parameter Real rhoc = 1.0e-7 "Central rest-mass density; small values are the Newtonian limit";
  parameter Boolean checkNewtonianLimit = true "Assert agreement with the Newtonian n = 1 profile (needs Gamma = 2 and rhoc*K << 1)";
  parameter Real tolerance = 1.0e-5 "Relative tolerance of the Lane-Emden comparison";
  parameter Real surfacePressure = 1.0e-12 "Pressure fraction pc*surfacePressure that defines the surface";
  final parameter Real n = 1.0 / (Gamma - 1.0) "Polytropic index";
  final parameter Real pc = GeneralRelativity.StressEnergy.polytropePressure(rhoc, K, Gamma) "Central pressure";
  final parameter Real ec = GeneralRelativity.StressEnergy.polytropeEnergyDensity(rhoc, K, Gamma) "Central energy density";
  final parameter Real alpha = sqrt((n + 1.0) * K * rhoc^(1.0 / n - 1.0) / (4.0 * pi)) "Lane-Emden length scale";
  final parameter Real massScale = 4.0 * pi * alpha^3 * rhoc "Lane-Emden mass scale";
  final parameter Real newtonianRadius = pi * alpha "n = 1 Lane-Emden radius";
  final parameter Real newtonianMass = pi * massScale "n = 1 Lane-Emden mass";
  Real r "Areal radius, equal to time";
  Real pHat "Pressure / central pressure";
  Real mHat "Gravitational mass inside r / massScale";
  Real baryonHat "Rest (baryon) mass inside r / massScale";
  Real pressure "Pressure";
  Real m "Gravitational mass function, g_rr = 1/(1 - 2m/r)";
  Real rho "Rest-mass density";
  Real energyDensity "Total energy density";
  Real xi "Lane-Emden radius r/alpha";
  Real thetaNewton "n = 1 Lane-Emden solution sin(xi)/xi";
  Real massNewton "Newtonian mass inside r";
  discrete Real radius(start = 0.0, fixed = true) "Surface radius, set when the pressure vanishes";
  discrete Real mass(start = 0.0, fixed = true) "Gravitational mass M, set at the surface";
  discrete Real baryonMass(start = 0.0, fixed = true) "Baryon mass, set at the surface";
initial equation
  // Series solution about the centre at the start radius r = StartTime (the equations are 0/0 at r = 0).
  pHat = 1.0 - 2.0 * pi / 3.0 * (ec + pc) * (ec + 3.0 * pc) * r^2 / pc;
  mHat = 4.0 * pi / 3.0 * ec * r^3 / massScale;
  baryonHat = 4.0 * pi / 3.0 * rhoc * r^3 / massScale;
equation
  r = time;
  pressure = pc * pHat;
  m = massScale * mHat;
  rho = GeneralRelativity.StressEnergy.polytropeDensity(pressure, K, Gamma);
  energyDensity = GeneralRelativity.StressEnergy.polytropeEnergyDensity(rho, K, Gamma);
  massScale * der(mHat) = 4.0 * pi * r^2 * energyDensity;
  massScale * der(baryonHat) = 4.0 * pi * r^2 * rho / sqrt(1.0 - 2.0 * m / r);
  pc * der(pHat) = -(energyDensity + pressure) * (m + 4.0 * pi * r^3 * pressure) / (r * (r - 2.0 * m));
  xi = r / alpha;
  thetaNewton = sin(xi) / xi;
  massNewton = massScale * (sin(xi) - xi * cos(xi));
  when pHat < surfacePressure then
    radius = r;
    mass = m;
    baryonMass = massScale * baryonHat;
    terminate("Reached the stellar surface");
  end when;
  if checkNewtonianLimit then
    assert(abs(pHat - max(thetaNewton, 0.0)^2) < tolerance or pHat < surfacePressure,
      "TOV pressure departs from the n = 1 Lane-Emden profile p = pc (sin(xi)/xi)^2");
    assert(abs(m - massNewton) < tolerance * newtonianMass or pHat < surfacePressure,
      "TOV mass departs from the n = 1 Lane-Emden mass 4 pi alpha^3 rhoc (sin(xi) - xi cos(xi))");
    when pHat < surfacePressure then
      assert(abs(r / newtonianRadius - 1.0) < tolerance, "Stellar radius is not the Lane-Emden radius pi alpha");
      assert(abs(m / newtonianMass - 1.0) < tolerance, "Stellar mass is not the Lane-Emden mass 4 pi^2 alpha^3 rhoc");
    end when;
  end if;
  assert(2.0 * m <= 8.0 / 9.0 * r, "Buchdahl bound 2M/R < 8/9 violated: no static star");
  annotation(experiment(StartTime = 4.0e-5, StopTime = 1.4, Tolerance = 1.0e-10, Interval = 1.0e-3),
    Documentation(info="<html>
<p>Static spherically symmetric perfect-fluid star with the polytropic equation of
state <code>p = K&rho;<sup>&Gamma;</sup></code> and energy density
<code>e = &rho; + p/(&Gamma; - 1)</code>
(<a href=\"modelica://GeneralRelativity.StressEnergy.polytropePressure\">StressEnergy.polytropePressure</a>,
<a href=\"modelica://GeneralRelativity.StressEnergy.polytropeEnergyDensity\">polytropeEnergyDensity</a>),
in units <code>G = c = 1</code>:</p>
<pre>dm/dr = 4&pi;r&sup2;e,   dp/dr = -(e + p)(m + 4&pi;r&sup3;p)/(r(r - 2m)),
dm<sub>b</sub>/dr = 4&pi;r&sup2;&rho;/sqrt(1 - 2m/r).</pre>
<p>The pressure and masses are integrated in units of the central pressure and of
the Lane-Emden mass scale <code>4&pi;&alpha;&sup3;&rho;<sub>c</sub></code>, with
<code>&alpha;&sup2; = (n + 1)K&rho;<sub>c</sub><sup>1/n - 1</sup>/(4&pi;)</code>, so the solver
tolerance is meaningful even when <code>p<sub>c</sub></code> is tiny. Integration starts
at a small radius <code>r<sub>0</sub></code> (the <code>StartTime</code>, 10<sup>-4</sup>&alpha; by default) from the
series solution <code>p = p<sub>c</sub> - (2&pi;/3)(e<sub>c</sub> + p<sub>c</sub>)(e<sub>c</sub> + 3p<sub>c</sub>)r<sub>0</sub>&sup2;</code>,
<code>m = (4&pi;/3)e<sub>c</sub>r<sub>0</sub>&sup3;</code>, and stops at
the surface where <code>p</code> drops below <code>surfacePressure&middot;p<sub>c</sub></code>;
<code>radius</code>, <code>mass</code> and <code>baryonMass</code> then hold the stellar
parameters. The exterior is
<a href=\"modelica://GeneralRelativity.Metrics.schwarzschild\">Schwarzschild</a> with
<code>r<sub>s</sub> = 2M</code>.</p>
<p><b>Check.</b> For <code>&Gamma; = 2</code> (<code>n = 1</code>) and
<code>K&rho;<sub>c</sub> &lt;&lt; 1</code> the star is Newtonian and the Lane-Emden solution is
<code>&theta; = sin&xi;/&xi;</code>: <code>p = p<sub>c</sub>&theta;&sup2;</code>,
<code>m = 4&pi;&alpha;&sup3;&rho;<sub>c</sub>(sin&xi; - &xi;cos&xi;)</code>,
<code>R = &pi;&alpha;</code>, <code>M = 4&pi;&sup2;&alpha;&sup3;&rho;<sub>c</sub></code>, independent of
<code>&rho;<sub>c</sub></code>. With the defaults (<code>K = 1</code>,
<code>&rho;<sub>c</sub> = 10<sup>-7</sup></code>, <code>2M/R &asymp; 4&times;10<sup>-7</sup></code>)
the relativistic corrections are below <code>10<sup>-6</sup></code>. The pressure profile is
asserted to <code>10<sup>-5</sup>p<sub>c</sub></code> and the mass profile to <code>10<sup>-5</sup>M</code>
throughout the run, and the radius and total mass to a relative <code>10<sup>-5</sup></code>
at the surface.</p>
<p><b>Relativistic star.</b> With <code>K = 100</code>, <code>&rho;<sub>c</sub> = 1.28&times;10<sup>-3</sup></code>,
<code>checkNewtonianLimit = false</code> and <code>StopTime = 12</code>, <code>StartTime = 4e-4</code> (lengths in units of
<code>GM<sub>sun</sub>/c&sup2;</code>) this is the standard numerical-relativity test star (Font et al. 2002); the model
gives <code>M = 1.4002</code>, baryon mass <code>1.5062</code> and areal radius
<code>R = 9.586</code>, so <code>2M/R = 0.29</code>.</p>
</html>"));
end PolytropeStar;
