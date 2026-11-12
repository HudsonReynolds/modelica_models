within GeneralRelativity.Perturbation;
function lenseThirringAcceleration "Gravitomagnetic (frame-dragging) acceleration of a test body about a spinning mass"
  input Real x[3](each unit = "m") "Position relative to the mass";
  input Real v[3](each unit = "m/s") "Coordinate velocity";
  input Real J[3](each unit = "kg.m2/s") "Spin angular momentum of the central body";
  input Real c(unit = "m/s") "Speed of light";
  input Real gravitationalConstant(unit = "m3/(kg.s2)") = GeneralRelativity.Constants.G "Newton's constant";
  output Real a[3](each unit = "m/s2") "Correction to add to the Newtonian acceleration";
protected
  Real r = sqrt(x * x);
algorithm
  a := 2.0 * gravitationalConstant / (c^2 * r^3) * (3.0 / r^2 * (x * J) * cross(x, v) + cross(v, J));
  annotation(Documentation(info="<html>
<p>IERS Conventions (2010) equation 10.12, Lense-Thirring term with
<code>&gamma; = 1</code>. It is the geodesic consequence of the
<code>g<sub>0i</sub></code> components in
<a href=\"modelica://GeneralRelativity.Metrics.earth\">Metrics.earth</a>: a
gyroscope or orbital node precesses about <code>J</code> at
<code>2GJ/(c&sup2;a&sup3;(1-e&sup2;)<sup>3/2</sup>)</code>, about 0.2 arcseconds per year
for LAGEOS.</p>
</html>"));
end lenseThirringAcceleration;
