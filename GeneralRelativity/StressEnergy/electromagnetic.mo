within GeneralRelativity.StressEnergy;
function electromagnetic "Maxwell stress-energy T_ab = (1/mu0) (F_ac F_b^c - 1/4 g_ab F_cd F^cd)"
  input Real g[4, 4] "Covariant metric";
  input Real F[4, 4] "Covariant field strength F_ab = d_a A_b - d_b A_a, antisymmetric";
  input Real mu0 = 4.0 * GeneralRelativity.Constants.pi "Vacuum permeability: 4 pi for geometrized Gaussian units, Modelica.Constants.mu_0 for SI";
  output Real T[4, 4] "T_ab";
protected
  Real gInv[4, 4];
  Boolean accepted;
  Real FMixed[4, 4] "F_a^c = F_ad g^dc";
  Real FUpper[4, 4] "F^ab";
  Real invariant "F_cd F^cd";
algorithm
  (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(g);
  FMixed := F * gInv;
  FUpper := gInv * FMixed;
  invariant := sum(F[a, b] * FUpper[a, b] for a in 1:4, b in 1:4);
  for a in 1:4 loop
    for b in 1:4 loop
      T[a, b] := (sum(F[a, k] * FMixed[b, k] for k in 1:4) - 0.25 * g[a, b] * invariant) / mu0;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>Units are fixed by <code>&mu;<sub>0</sub></code>:</p>
<ul>
<li><b>Geometrized Gaussian</b> (default, <code>G = c = 4&pi;&epsilon;<sub>0</sub> = 1</code>,
<code>&mu;<sub>0</sub> = 4&pi;</code>):
<code>T<sub>ab</sub> = (1/4&pi;)(F<sub>ac</sub>F<sub>b</sub><sup>c</sup> - &frac14;g<sub>ab</sub>F<sub>cd</sub>F<sup>cd</sup>)</code>.
A point charge has <code>F<sub>tr</sub> = Q/r&sup2;</code> and the field equations are
<code>G<sub>ab</sub> = 8&pi;T<sub>ab</sub></code>, which is solved by
<a href=\"modelica://GeneralRelativity.Metrics.reissnerNordstrom\">Metrics.reissnerNordstrom</a>
with <code>r<sub>Q</sub> = Q</code>.</li>
<li><b>SI</b>: pass <code>&mu;<sub>0</sub> = 4&pi;&times;10<sup>-7</sup></code> N/A&sup2; and
<code>F<sub>0i</sub> = -E<sub>i</sub>/c</code>, <code>F<sub>ij</sub> = &epsilon;<sub>ijk</sub>B<sub>k</sub></code>
(<a href=\"modelica://GeneralRelativity.StressEnergy.fieldStrength\">fieldStrength</a>).
Then <code>T<sub>00</sub> = &epsilon;<sub>0</sub>E&sup2;/2 + B&sup2;/(2&mu;<sub>0</sub>)</code> in
J/m&sup3; and <code>T<sub>0i</sub> = -S<sub>i</sub>/c</code> with the Poynting vector
<code>S = E &times; B/&mu;<sub>0</sub></code>; use the field equations with
<code>8&pi;G/c<sup>4</sup></code>.</li>
</ul>
<p>The tensor is traceless in four dimensions, so electrovacuum metrics have
<code>R = 0</code>, and it satisfies the dominant energy condition.</p>
</html>"));
end electromagnetic;
