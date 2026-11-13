within GeneralRelativity.Metrics;
function kerrSchild "Kerr-Newman family in Cartesian Kerr-Schild coordinates {ct, x, y, z}, regular across the horizon; p = {rs}, {rs, a} or {rs, a, rQ}"
  extends GeneralRelativity.Interfaces.MetricFunction;
protected
  Real rs = p[1] "Schwarzschild radius 2GM/c^2";
  Real a = if size(p, 1) >= 2 then p[2] else 0.0 "Spin length J/(Mc), about +z";
  Real rQ = if size(p, 1) >= 3 then p[3] else 0.0 "Charge length, rQ^2 = G Q^2/(4 pi eps0 c^4)";
  Real rho2 = x[2]^2 + x[3]^2 + x[4]^2;
  Real w = rho2 - a^2;
  Real r = sqrt(0.5 * (w + sqrt(w^2 + 4.0 * a^2 * x[4]^2))) "Boyer-Lindquist radius";
  Real f = r^2 * (rs * r - rQ^2) / (r^4 + a^2 * x[4]^2);
  Real k[4] = {1.0, (r * x[2] + a * x[3]) / (r^2 + a^2), (r * x[3] - a * x[2]) / (r^2 + a^2), x[4] / r};
algorithm
  g := diagonal({-1.0, 1.0, 1.0, 1.0});
  for i in 1:4 loop
    for j in 1:4 loop
      g[i, j] := g[i, j] + f * k[i] * k[j];
    end for;
  end for;
  annotation(Documentation(info="<html>
<p><code>g<sub>&mu;&nu;</sub> = &eta;<sub>&mu;&nu;</sub> + f k<sub>&mu;</sub>k<sub>&nu;</sub></code> with the
null covector <code>k = (1, (rx + ay)/(r&sup2; + a&sup2;), (ry &minus; ax)/(r&sup2; + a&sup2;), z/r)</code> and
<code>f = r&sup2;(r<sub>s</sub>r &minus; r<sub>Q</sub>&sup2;)/(r&sup4; + a&sup2;z&sup2;)</code>, where the Boyer-Lindquist
radius solves <code>r&sup4; &minus; (x&sup2; + y&sup2; + z&sup2; &minus; a&sup2;)r&sup2; &minus; a&sup2;z&sup2; = 0</code>. Since
<code>k</code> is null for <code>&eta;</code>, the inverse is exactly
<code>g<sup>&mu;&nu;</sup> = &eta;<sup>&mu;&nu;</sup> &minus; f k<sup>&mu;</sup>k<sup>&nu;</sup></code>.</p>
<p>Special cases: <code>a = r<sub>Q</sub> = 0</code> is Schwarzschild in ingoing Eddington-Finkelstein form,
<code>r<sub>Q</sub> = 0</code> Kerr, <code>a = 0</code> Reissner-Nordstr&ouml;m. The sense of rotation and the
curvature invariants agree with <a href=\"modelica://GeneralRelativity.Metrics.kerr\">kerr</a>, as
<code>Tests.GeneralRelativityVisualizationTests</code> checks. The Three.js ray tracer integrates null
geodesics of exactly this metric, because it has no coordinate singularity at the horizon. The
Boyer-Lindquist point <code>(r, &theta;, &phi;)</code> sits at <code>x + iy = (r + ia) sin&theta; e<sup>i&phi;&#771;</sup></code>,
<code>z = r cos&theta;</code>, where <code>&phi;&#771;</code> differs from <code>&phi;</code> by a function of <code>r</code>.</p>
</html>"));
end kerrSchild;
