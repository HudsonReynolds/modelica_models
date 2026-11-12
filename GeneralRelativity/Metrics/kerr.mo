within GeneralRelativity.Metrics;
function kerr "Kerr metric in Boyer-Lindquist coordinates {ct, r, theta, phi}; p = {rs, a}"
  extends GeneralRelativity.Interfaces.MetricFunction;
protected
  Real rs = p[1] "Schwarzschild radius 2GM/c^2";
  Real a = p[2] "Spin parameter J/(M c), a length";
  Real r = x[2];
  Real sinTheta = sin(x[3]);
  Real cosTheta = cos(x[3]);
  Real Sigma = r^2 + a^2 * cosTheta^2;
  Real Delta = r^2 - rs * r + a^2;
algorithm
  g := zeros(4, 4);
  g[1, 1] := -(1.0 - rs * r / Sigma);
  g[1, 4] := -rs * r * a * sinTheta^2 / Sigma;
  g[4, 1] := g[1, 4];
  g[2, 2] := Sigma / Delta;
  g[3, 3] := Sigma;
  g[4, 4] := (r^2 + a^2 + rs * r * a^2 * sinTheta^2 / Sigma) * sinTheta^2;
  annotation(Documentation(info="<html>
<p>Rotating vacuum black hole with <code>a = J/(Mc)</code>; <code>a = 0</code>
recovers Schwarzschild. The outer horizon is at
<code>r<sub>+</sub> = r<sub>s</sub>/2 + sqrt(r<sub>s</sub>&sup2;/4 - a&sup2;)</code>, so
<code>|a| &le; r<sub>s</sub>/2</code>. The <code>g<sub>t&phi;</sub></code> term is
frame dragging: zero-angular-momentum observers rotate with
<code>&omega; = -g<sub>t&phi;</sub>/g<sub>&phi;&phi;</sub></code>.</p>
</html>"));
end kerr;
