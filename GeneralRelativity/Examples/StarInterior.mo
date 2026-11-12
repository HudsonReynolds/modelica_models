within GeneralRelativity.Examples;
model StarInterior "Tolman-Oppenheimer-Volkoff integration of a constant-density star against the interior Schwarzschild solution; Modelica time is the radius"
  parameter Real rho = 0.05 "Energy density in geometric units; must satisfy 2M/R < 8/9";
  parameter Real R = 1.0 "Stellar radius";
  parameter Real r0 = 1.0e-3 "Start radius; the equations are regular but 0/0 at the centre";
  final parameter Real M = 4.0 / 3.0 * GeneralRelativity.Constants.pi * rho * R^3 "Gravitational mass";
  final parameter Real compactness = 2.0 * M / R;
  Real r "Radius, equal to time";
  Real m(start = 4.0 / 3.0 * GeneralRelativity.Constants.pi * rho * r0^3, fixed = true) "Mass inside r";
  Real pressure(start = exactPressure(r0, rho, M, R), fixed = true) "Pressure from the TOV equation";
  Real nu(start = 0.0, fixed = true) "Metric potential: g_tt = -exp(2 nu), up to a constant fixed at the surface";
  Real pressureExact "Interior Schwarzschild pressure";
  Real g[4, 4] "Interior metric diag(-exp(2 nu), 1/(1 - 2m/r), r^2, r^2) at the equator";
  Real u[4] "Static fluid four-velocity";
  Real T[4, 4] "Perfect-fluid stress-energy T_ab";
  Real energyDensityCheck "T_ab u^a u^b, equal to rho";
equation
  r = time;
  der(m) = 4.0 * GeneralRelativity.Constants.pi * r^2 * rho;
  der(pressure) = -(rho + pressure) * (m + 4.0 * GeneralRelativity.Constants.pi * r^3 * pressure) / (r * (r - 2.0 * m));
  der(nu) = (m + 4.0 * GeneralRelativity.Constants.pi * r^3 * pressure) / (r * (r - 2.0 * m));
  pressureExact = exactPressure(r, rho, M, R);
  g = diagonal({-exp(2.0 * nu), 1.0 / (1.0 - 2.0 * m / r), r^2, r^2});
  u = {exp(-nu), 0.0, 0.0, 0.0};
  T = GeneralRelativity.StressEnergy.perfectFluid(g, u, rho, pressure, 1.0);
  energyDensityCheck = GeneralRelativity.StressEnergy.energyDensity(T, u, 1.0);
  assert(compactness < 8.0 / 9.0, "Buchdahl limit violated: no static solution exists");
  assert(abs(pressure - pressureExact) < 1.0e-6 * rho, "TOV pressure disagrees with the interior Schwarzschild solution");
  assert(abs(energyDensityCheck - rho) < 1.0e-9 * rho, "Perfect-fluid energy density is not rho in the fluid frame");
protected
  function exactPressure "Interior Schwarzschild solution for constant density"
    input Real r;
    input Real rho;
    input Real M;
    input Real R;
    output Real p;
  protected
    Real surface = sqrt(1.0 - 2.0 * M / R);
    Real local = sqrt(1.0 - 2.0 * M * r^2 / R^3);
  algorithm
    p := rho * (local - surface) / (3.0 * surface - local);
  end exactPressure;
  annotation(experiment(StartTime = 1.0e-3, StopTime = 0.999, Tolerance = 1.0e-10, Interval = 1.0e-3),
    Documentation(info="<html>
<p>The static spherically symmetric field equations with a perfect fluid reduce
to <code>dm/dr = 4&pi;r&sup2;&rho;</code> and the TOV equation
<code>dp/dr = -(&rho; + p)(m + 4&pi;r&sup3;p)/(r(r - 2m))</code> in units
<code>G = c = 1</code>. For constant density the pressure is known in closed form
and is reproduced to the solver tolerance. The compactness must stay below the
Buchdahl bound <code>8/9</code>. Start the simulation at <code>r<sub>0</sub></code>;
the model's <code>experiment</code> annotation sets <code>StartTime = 1e-3</code>.</p>
</html>"));
end StarInterior;
