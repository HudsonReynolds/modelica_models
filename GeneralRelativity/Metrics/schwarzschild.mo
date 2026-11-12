within GeneralRelativity.Metrics;
function schwarzschild "Schwarzschild metric in {ct, r, theta, phi}; p = {rs}"
  extends GeneralRelativity.Interfaces.MetricFunction;
protected
  Real rs = p[1] "Schwarzschild radius 2GM/c^2";
  Real r = x[2];
  Real f = 1.0 - rs / r;
algorithm
  g := zeros(4, 4);
  g[1, 1] := -f;
  g[2, 2] := 1.0 / f;
  g[3, 3] := r^2;
  g[4, 4] := r^2 * sin(x[3])^2;
  annotation(Documentation(info="<html>
<p><code>ds&sup2; = -(1 - r<sub>s</sub>/r) c&sup2;dt&sup2; + dr&sup2;/(1 - r<sub>s</sub>/r)
+ r&sup2;(d&theta;&sup2; + sin&sup2;&theta; d&phi;&sup2;)</code>. Vacuum: Ricci vanishes and
the Kretschmann scalar is <code>12 r<sub>s</sub>&sup2;/r<sup>6</sup></code>. The chart is
singular at <code>r = r<sub>s</sub></code>; the exterior is <code>r &gt; r<sub>s</sub></code>.</p>
<p>For a static observer the conserved specific energy is
<code>E = -g<sub>tt</sub> u<sup>t</sup></code> and the specific angular momentum is
<code>L = g<sub>&phi;&phi;</sub> u<sup>&phi;</sup></code>; see
<a href=\"modelica://GeneralRelativity.Geodesics.killingConserved\">Geodesics.killingConserved</a>.</p>
</html>"));
end schwarzschild;
