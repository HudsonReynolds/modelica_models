within GeneralRelativity.Metrics;
function reissnerNordstrom "Reissner-Nordstrom charged black hole in {ct, r, theta, phi}; p = {rs, rQ}"
  extends GeneralRelativity.Interfaces.MetricFunction;
protected
  Real rs = p[1] "Schwarzschild radius 2GM/c^2";
  Real rQ = p[2] "Charge length, rQ^2 = G Q^2 / (4 pi epsilon0 c^4)";
  Real r = x[2];
  Real f = 1.0 - rs / r + rQ^2 / r^2;
algorithm
  g := zeros(4, 4);
  g[1, 1] := -f;
  g[2, 2] := 1.0 / f;
  g[3, 3] := r^2;
  g[4, 4] := r^2 * sin(x[3])^2;
  annotation(Documentation(info="<html>
<p><code>ds&sup2; = -f c&sup2;dt&sup2; + dr&sup2;/f + r&sup2;(d&theta;&sup2; + sin&sup2;&theta; d&phi;&sup2;)</code>
with <code>f = 1 - r<sub>s</sub>/r + r<sub>Q</sub>&sup2;/r&sup2;</code>.
The charge enters as the length
<code>r<sub>Q</sub>&sup2; = G Q&sup2;/(4&pi;&epsilon;<sub>0</sub>c<sup>4</sup>)</code> (SI); in
geometrized Gaussian units (<code>G = c = 4&pi;&epsilon;<sub>0</sub> = 1</code>)
<code>r<sub>Q</sub> = Q</code> and <code>r<sub>s</sub> = 2M</code>.</p>
<p>The metric solves the Einstein-Maxwell equations with the radial Coulomb field
<code>F<sub>tr</sub> = Q/r&sup2;</code> (geometrized Gaussian): its Einstein tensor equals
<code>8&pi;T<sub>ab</sub></code> of
<a href=\"modelica://GeneralRelativity.StressEnergy.electromagnetic\">StressEnergy.electromagnetic</a>
and the Ricci scalar vanishes because that source is traceless. Horizons sit at
<code>r<sub>&plusmn;</sub> = r<sub>s</sub>/2 &plusmn; sqrt(r<sub>s</sub>&sup2;/4 - r<sub>Q</sub>&sup2;)</code>;
<code>r<sub>Q</sub> = 0</code> recovers
<a href=\"modelica://GeneralRelativity.Metrics.schwarzschild\">schwarzschild</a>.</p>
</html>"));
end reissnerNordstrom;
