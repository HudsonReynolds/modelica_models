within GeneralRelativity.NewtonCartan;
function pointMassSpherical "Newton-Cartan connection of a point mass in {ct, r, theta, phi}; p = {rs}"
  extends GeneralRelativity.Interfaces.ConnectionFunction;
protected
  Real rs = p[1] "2GM/c^2, so that GM/(c^2 r^2) = rs/(2 r^2)";
  Real r = x[2];
  Real theta = x[3];
algorithm
  Gamma := zeros(4, 4, 4);
  // Gravitational part: Gamma^r_00 = (dPhi/dr)/c^2 with Phi = -GM/r.
  Gamma[2, 1, 1] := rs / (2.0 * r^2);
  // Flat Euclidean connection of the spherical chart.
  Gamma[2, 3, 3] := -r;
  Gamma[2, 4, 4] := -r * sin(theta)^2;
  Gamma[3, 2, 3] := 1.0 / r;
  Gamma[3, 3, 2] := 1.0 / r;
  Gamma[3, 4, 4] := -sin(theta) * cos(theta);
  Gamma[4, 2, 4] := 1.0 / r;
  Gamma[4, 4, 2] := 1.0 / r;
  Gamma[4, 3, 4] := cos(theta) / sin(theta);
  Gamma[4, 4, 3] := Gamma[4, 3, 4];
  annotation(Documentation(info="<html>
<p>The chart and parameter layout match
<a href=\"modelica://GeneralRelativity.Metrics.schwarzschild\">Metrics.schwarzschild</a>, so a
Newtonian body can be integrated beside a Schwarzschild geodesic from the same
initial data. Compared with
<a href=\"modelica://GeneralRelativity.Metrics.schwarzschildChristoffel\">Metrics.schwarzschildChristoffel</a>,
<code>&Gamma;<sup>r</sup><sub>00</sub></code> lacks the factor <code>1 &minus; r<sub>s</sub>/r</code>
and the <code>r<sub>s</sub></code>-dependent spatial and <code>&Gamma;<sup>0</sup><sub>0r</sub></code>
terms are absent: absolute time does not depend on position.</p>
</html>"));
end pointMassSpherical;
