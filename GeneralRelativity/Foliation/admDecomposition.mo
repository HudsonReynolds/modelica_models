within GeneralRelativity.Foliation;
function admDecomposition "Lapse, shift and spatial metric of a covariant 4-metric"
  input Real g[4, 4] "Covariant metric; x[1] = c*t";
  output Real alpha "Lapse";
  output Real betaUp[3] "Shift beta^i";
  output Real betaDown[3] "Shift beta_i = g_0i";
  output Real gamma[3, 3] "Spatial metric gamma_ij";
  output Real gammaInv[3, 3] "Inverse spatial metric gamma^ij";
  output Boolean accepted "false when gamma is singular or the slice x^0 = const is not spacelike";
protected
  Real alpha2;
  Boolean inverted;
algorithm
  gamma := g[2:4, 2:4];
  (gammaInv, inverted) := GeneralRelativity.Tensors.inverseMetric(gamma);
  betaDown := g[1, 2:4];
  betaUp := gammaInv * betaDown;
  alpha2 := betaUp * betaDown - g[1, 1];
  accepted := inverted and alpha2 > 0.0;
  alpha := sqrt(max(alpha2, 0.0));
  annotation(Documentation(info="<html>
<p>Pure algebra on <code>g<sub>00</sub> = -&alpha;&sup2; + &beta;<sub>k</sub>&beta;<sup>k</sup></code>,
<code>g<sub>0i</sub> = &beta;<sub>i</sub></code>, <code>g<sub>ij</sub> = &gamma;<sub>ij</sub></code>.
Equivalently <code>&alpha; = (-g<sup>00</sup>)<sup>-1/2</sup></code> and
<code>&beta;<sup>i</sup> = -g<sup>0i</sup>/g<sup>00</sup></code>. A non-positive
<code>&alpha;&sup2;</code> means the surface <code>x<sup>0</sup> = const</code> is not
spacelike, for example inside a Schwarzschild horizon in Schwarzschild time.</p>
</html>"));
end admDecomposition;
