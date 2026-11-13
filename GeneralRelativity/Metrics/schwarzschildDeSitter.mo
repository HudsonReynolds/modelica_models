within GeneralRelativity.Metrics;
function schwarzschildDeSitter "Schwarzschild-de Sitter (Kottler) metric in {ct, r, theta, phi}; p = {rs, Lambda}"
  extends GeneralRelativity.Interfaces.MetricFunction;
protected
  Real rs = p[1] "Schwarzschild radius 2GM/c^2";
  Real Lambda = p[2] "Cosmological constant (1/length^2)";
  Real r = x[2];
  Real f = 1.0 - rs / r - Lambda * r^2 / 3.0;
algorithm
  g := zeros(4, 4);
  g[1, 1] := -f;
  g[2, 2] := 1.0 / f;
  g[3, 3] := r^2;
  g[4, 4] := r^2 * sin(x[3])^2;
  annotation(Documentation(info="<html>
<p>Static chart with <code>f = 1 - r<sub>s</sub>/r - &Lambda;r&sup2;/3</code>. It solves
<code>G<sub>ab</sub> + &Lambda;g<sub>ab</sub> = 0</code>, equivalently
<code>R<sub>ab</sub> = &Lambda;g<sub>ab</sub></code> and <code>R = 4&Lambda;</code>.
<code>r<sub>s</sub> = 0</code> is static de Sitter space with the cosmological horizon at
<code>r = sqrt(3/&Lambda;)</code>; <code>&Lambda; &lt; 0</code> gives
Schwarzschild-anti-de Sitter. The static region lies between the black-hole and
cosmological horizons, the two positive roots of <code>f</code>.</p>
</html>"));
end schwarzschildDeSitter;
