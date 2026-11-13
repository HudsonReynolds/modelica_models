within GeneralRelativity.Foliation;
function normalVector "Future-pointing unit normal of the slices x^0 = const"
  input Real g[4, 4] "Covariant metric; x[1] = c*t";
  output Real nUp[4] "n^a = (1, -beta^i)/alpha";
  output Real nDown[4] "n_a = (-alpha, 0, 0, 0)";
protected
  Real alpha;
  Real betaUp[3];
  Real betaDown[3];
  Real gamma[3, 3];
  Real gammaInv[3, 3];
  Boolean accepted;
algorithm
  (alpha, betaUp, betaDown, gamma, gammaInv, accepted) := admDecomposition(g);
  nDown := {-alpha, 0.0, 0.0, 0.0};
  nUp := cat(1, {1.0}, -betaUp) / alpha;
  annotation(Documentation(info="<html>
<p><code>g<sub>ab</sub>n<sup>a</sup>n<sup>b</sup> = -1</code> and
<code>n<sub>a</sub> = g<sub>ab</sub>n<sup>b</sup></code>; the normal points to
increasing <code>x<sup>0</sup></code>. These are the four-velocities of the
Eulerian (normal) observers.</p>
</html>"));
end normalVector;
