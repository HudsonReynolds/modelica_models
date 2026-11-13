within GeneralRelativity.Harmonics;
function gaussLegendre "Gauss-Legendre nodes and weights on [-1, 1]"
  input Integer n "Number of nodes, n >= 1";
  output Real x[n] "Nodes in ascending order";
  output Real w[n] "Weights, summing to 2";
protected
  Real z;
  Real zOld;
  Real p0;
  Real p1;
  Real p2;
  Real dp "P_n'(z)";
  Integer iterations;
algorithm
  for i in 1:n loop
    // Initial guess near the i-th largest root, then Newton on P_n.
    z := cos(GeneralRelativity.Constants.pi * (i - 0.25) / (n + 0.5));
    zOld := 2.0;
    iterations := 0;
    p1 := 1.0;
    dp := 1.0;
    while abs(z - zOld) > 1.0e-15 and iterations < 100 loop
      p0 := 0.0;
      p1 := 1.0;
      for k in 1:n loop
        p2 := p0;
        p0 := p1;
        p1 := ((2 * k - 1) * z * p0 - (k - 1) * p2) / k;
      end for;
      // p1 = P_n(z), p0 = P_(n-1)(z)
      dp := n * (z * p1 - p0) / (z^2 - 1.0);
      zOld := z;
      z := zOld - p1 / dp;
      iterations := iterations + 1;
    end while;
    x[n + 1 - i] := z;
    w[n + 1 - i] := 2.0 / ((1.0 - z^2) * dp^2);
  end for;
  annotation(Documentation(info="<html>
<p>The nodes are the roots of the Legendre polynomial <code>P<sub>n</sub></code>, found by
Newton's method from the asymptotic guess <code>cos(&pi;(i - &frac14;)/(n + &frac12;))</code>; the
weights are <code>2/((1 - x&sup2;)P<sub>n</sub>'(x)&sup2;)</code>. The rule integrates
polynomials of degree <code>2n - 1</code> exactly. On the sphere, nodes in
<code>cos&theta;</code> combined with <code>N<sub>&phi;</sub></code> equally spaced azimuths
(weight <code>2&pi;/N<sub>&phi;</sub></code>) integrate every product of harmonics of total
degree below <code>2n</code> and total order below <code>N<sub>&phi;</sub></code> exactly.</p>
</html>"));
end gaussLegendre;
