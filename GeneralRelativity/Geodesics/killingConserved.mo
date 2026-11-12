within GeneralRelativity.Geodesics;
function killingConserved "Constant of motion g_ab xi^a u^b along a geodesic for a Killing vector xi"
  input Real g[4, 4] "Covariant metric";
  input Real xi[4] "Killing vector components; {1,0,0,0} for energy, {0,0,0,1} for axial angular momentum";
  input Real u[4] "Four-velocity";
  output Real constant_ "g(xi, u)";
algorithm
  constant_ := xi * (g * u);
  annotation(Inline = true, Documentation(info="<html>
<p>In Schwarzschild coordinates with <code>x[1] = ct</code>, the specific energy
per <code>c</code> is <code>E = -g(&part;<sub>t</sub>, u) = (1 - r<sub>s</sub>/r) u<sup>t</sup></code>
and the specific angular momentum is <code>L = g(&part;<sub>&phi;</sub>, u) = r&sup2;sin&sup2;&theta; u<sup>&phi;</sup></code>.</p>
</html>"));
end killingConserved;
