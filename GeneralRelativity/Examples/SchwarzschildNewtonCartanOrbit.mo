within GeneralRelativity.Examples;
model SchwarzschildNewtonCartanOrbit "Schwarzschild geodesic beside the Newton-Cartan geodesic from the same initial data, advanced in the same coordinate time"
  extends GeneralRelativity.Examples.SchwarzschildOrbit;
  final parameter Real omega0 = u_start[4] / u_start[1] "Initial coordinate angular velocity d phi / d(ct), shared by both bodies";
  GeneralRelativity.Geodesics.AffineGeodesic newton(
    redeclare function connection = GeneralRelativity.NewtonCartan.pointMassSpherical,
    p = {rs},
    c = c,
    x_start = x_start,
    u_start = {1.0, 0.0, 0.0, omega0},
    lambdaRate = u[1]) "Newton-Cartan body; its affine parameter is ct, advanced at c dt/d tau of the relativistic body";
  final parameter Real newtonEnergyExpected = 0.5 * (r0 * omega0)^2 - rs / (2.0 * r0)
    "Newtonian specific energy over c^2 at the start";
  final parameter Real newtonAngularMomentumExpected = r0^2 * omega0 "Newtonian specific angular momentum over c";
  final parameter Real newtonSemiMajorAxis = -rs / (4.0 * newtonEnergyExpected) "a = -GM / (2 E)";
  final parameter Real newtonPeriodExpected = 2.0 * GeneralRelativity.Constants.pi * sqrt(newtonSemiMajorAxis^3 / M)
    "Kepler period in units of ct: 2 pi sqrt(a^3 / M) with M = GM/c^2";
  Real newtonEnergy "v^2/2 - GM/r over c^2; constant along the Newton-Cartan geodesic";
  Real newtonAngularMomentum "r^2 d phi / d(ct); constant";
  Real separation "Distance between the two bodies in the equatorial (r, phi) plane drawn as flat";
  Real coordinateTimeLag "x[1] of the Newtonian body minus x[1] of the relativistic body; zero by construction";
  discrete Real newtonPeriapsisPhi(start = 0.0, fixed = true) "Azimuth of the last Newtonian periapsis";
  discrete Real newtonPeriapsisTime(start = 0.0, fixed = true) "ct at the last Newtonian periapsis";
  discrete Real newtonPeriapsisAdvance(start = 0.0, fixed = true) "Azimuth gained per Newtonian orbit beyond 2 pi; zero for a closed Kepler ellipse";
  discrete Real newtonPeriod(start = 0.0, fixed = true) "ct between Newtonian periapsides";
  discrete Integer newtonPeriapsides(start = 0, fixed = true) "Newtonian periapsis passages seen";
equation
  newtonEnergy = 0.5 * (newton.u[2]^2 + (newton.x[2] * newton.u[3])^2 + (newton.x[2] * sin(newton.x[3]) * newton.u[4])^2)
    - rs / (2.0 * newton.x[2]);
  newtonAngularMomentum = (newton.x[2] * sin(newton.x[3]))^2 * newton.u[4];
  separation = sqrt((x[2] * cos(x[4]) - newton.x[2] * cos(newton.x[4]))^2 + (x[2] * sin(x[4]) - newton.x[2] * sin(newton.x[4]))^2);
  coordinateTimeLag = newton.x[1] - x[1];
  when newton.u[2] > 0.0 then
    newtonPeriapsisPhi = newton.x[4];
    newtonPeriapsisTime = newton.x[1];
    newtonPeriapsides = pre(newtonPeriapsides) + 1;
    newtonPeriapsisAdvance = if pre(newtonPeriapsides) >= 1 then newton.x[4] - pre(newtonPeriapsisPhi) - 2.0 * GeneralRelativity.Constants.pi else 0.0;
    newtonPeriod = if pre(newtonPeriapsides) >= 1 then newton.x[1] - pre(newtonPeriapsisTime) else 0.0;
  end when;
  assert(abs(newtonEnergy - newtonEnergyExpected) < 1.0e-6 * abs(newtonEnergyExpected), "Newton-Cartan geodesic lost Newtonian energy");
  assert(abs(newtonAngularMomentum - newtonAngularMomentumExpected) < 1.0e-6 * newtonAngularMomentumExpected,
    "Newton-Cartan geodesic lost angular momentum");
  assert(abs(coordinateTimeLag) < 1.0e-6 * (1.0 + abs(x[1])), "The two bodies are no longer at the same coordinate time");
  assert(abs(newtonPeriapsisAdvance) < 1.0e-5, "The Newton-Cartan orbit precesses; a Kepler ellipse must close");
  assert(newtonPeriapsides < 2 or abs(newtonPeriod / newtonPeriodExpected - 1.0) < 1.0e-6,
    "Newton-Cartan orbital period is not Kepler's 2 pi sqrt(a^3 / GM)");
  annotation(experiment(StartTime = 0.0, StopTime = 1500.0, Tolerance = 1.0e-10, Interval = 0.5),
    Documentation(info="<html>
<p>Two free-falling bodies leave periastron <code>r<sub>0</sub> = 10M</code> together:</p>
<ul>
<li>the relativistic body follows a geodesic of the
<a href=\"modelica://GeneralRelativity.Metrics.schwarzschild\">Schwarzschild</a> metric
(everything inherited from
<a href=\"modelica://GeneralRelativity.Examples.SchwarzschildOrbit\">SchwarzschildOrbit</a>, including its checks);</li>
<li>the Newtonian body <code>newton</code> follows a geodesic of the
<a href=\"modelica://GeneralRelativity.NewtonCartan.pointMassSpherical\">Newton-Cartan connection</a>
of the same mass in the same chart.</li>
</ul>
<p>Both start with the same position and the same coordinate velocity
<code>dx/dt</code>, where <code>t</code> is Schwarzschild coordinate time, the time of static
observers far away, which Newton-Cartan theory takes as absolute time. Modelica
<code>time</code> is the proper time of the relativistic body; the Newtonian body uses
<code>&lambda; = ct</code> and is advanced at <code>d&lambda;/d&tau; = u<sup>t</sup></code>, so at every
instant the two bodies sit at the same coordinate time
(<code>coordinateTimeLag</code> stays zero) and their positions can be compared directly.</p>
<p>What the comparison shows, with the default <code>L = 4.2M</code>:</p>
<ul>
<li>the Newtonian orbit is a closed Kepler ellipse: <code>newtonPeriapsisAdvance</code>
stays zero and the period equals <code>2&pi;&radic;(a&sup3;/GM)</code>, both asserted;</li>
<li>the relativistic orbit precesses by about 2.1 rad per orbit and reaches a larger
apoapsis (about <code>21M</code> against about <code>15M</code>): at the same coordinate
velocity gravity is stronger near the hole than Newton's law, and
<code>&Gamma;<sup>r</sup><sub>00</sub></code> differs by the factor <code>1 &minus; r<sub>s</sub>/r</code>
while <code>&Gamma;<sup>0</sup><sub>0r</sub></code> adds gravitational time dilation;</li>
<li><code>separation</code> grows from zero as the two pictures of gravity disagree.</li>
</ul>
<p><a href=\"modelica://GeneralRelativity.Examples.SchwarzschildOrbitScene\">SchwarzschildOrbitScene</a>
writes both trajectories to the Three.js viewer, which animates them together.</p>
</html>"));
end SchwarzschildNewtonCartanOrbit;
