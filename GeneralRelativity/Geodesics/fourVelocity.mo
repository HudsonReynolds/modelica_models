within GeneralRelativity.Geodesics;
function fourVelocity "Normalized four-velocity from a coordinate velocity dx^i/dt: g(u,u) = -c^2 for massive particles"
  input Real g[4, 4] "Covariant metric";
  input Real v[3] "Spatial coordinate velocity dx^i/dt in chart units per second";
  input Real c(unit = "m/s") = 1.0 "Speed of light";
  output Real u[4] "Contravariant four-velocity dx/d tau with x[1] = c t";
protected
  Real w[4] = {1.0, v[1] / c, v[2] / c, v[3] / c} "dx/d(ct)";
  Real norm = -(w * (g * w)) "-(g(w,w)); positive inside the light cone";
algorithm
  assert(norm > 0.0, "fourVelocity: the coordinate velocity is not timelike");
  u := c * w / sqrt(norm);
  annotation(Documentation(info="<html>
<p><code>u<sup>0</sup> = c dt/d&tau;</code> is the chart Lorentz factor times
<code>c</code>; in Minkowski space <code>u = &gamma;(c, v)</code>. For null curves
build <code>u</code> directly with
<a href=\"modelica://GeneralRelativity.Geodesics.normResidual\">normResidual</a>
as the check.</p>
</html>"));
end fourVelocity;
