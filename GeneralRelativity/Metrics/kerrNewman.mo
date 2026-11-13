within GeneralRelativity.Metrics;
function kerrNewman "Kerr-Newman charged rotating black hole in Boyer-Lindquist coordinates {ct, r, theta, phi}; p = {rs, a, rQ}"
  extends GeneralRelativity.Interfaces.MetricFunction;
protected
  Real rs = p[1] "Schwarzschild radius 2GM/c^2";
  Real a = p[2] "Spin parameter J/(M c), a length";
  Real rQ = p[3] "Charge length, rQ^2 = G Q^2 / (4 pi epsilon0 c^4)";
  Real r = x[2];
  Real sinTheta = sin(x[3]);
  Real cosTheta = cos(x[3]);
  Real Sigma = r^2 + a^2 * cosTheta^2;
  Real Delta = r^2 - rs * r + a^2 + rQ^2;
  Real w = rs * r - rQ^2 "Mass-minus-charge combination that replaces rs*r of Kerr";
algorithm
  g := zeros(4, 4);
  g[1, 1] := -(1.0 - w / Sigma);
  g[1, 4] := -w * a * sinTheta^2 / Sigma;
  g[4, 1] := g[1, 4];
  g[2, 2] := Sigma / Delta;
  g[3, 3] := Sigma;
  g[4, 4] := (r^2 + a^2 + w * a^2 * sinTheta^2 / Sigma) * sinTheta^2;
  annotation(Documentation(info="<html>
<p>Kerr with <code>r<sub>s</sub> r</code> replaced by
<code>r<sub>s</sub> r - r<sub>Q</sub>&sup2;</code>, so
<code>&Delta; = r&sup2; - r<sub>s</sub> r + a&sup2; + r<sub>Q</sub>&sup2;</code> and
<code>&Sigma; = r&sup2; + a&sup2;cos&sup2;&theta;</code>. The charge length is defined as in
<a href=\"modelica://GeneralRelativity.Metrics.reissnerNordstrom\">reissnerNordstrom</a>.
<code>a = 0</code> gives Reissner-Nordstrom and <code>r<sub>Q</sub> = 0</code> gives
<a href=\"modelica://GeneralRelativity.Metrics.kerr\">kerr</a>.</p>
<p>The source is the field of the vector potential (geometrized Gaussian units)
<code>A = -(Q r/&Sigma;)(c dt - a sin&sup2;&theta; d&phi;)</code>, i.e.
<code>F = dA</code>; the Ricci scalar vanishes. Horizons exist for
<code>a&sup2; + r<sub>Q</sub>&sup2; &le; r<sub>s</sub>&sup2;/4</code>.</p>
</html>"));
end kerrNewman;
