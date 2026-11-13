within GeneralRelativity.Harmonics;
function gaunt "Gaunt coefficient: integral of Y_l1m1 Y_l2m2 Y_l3m3 over the unit sphere"
  input Integer l1 "Degree of the first harmonic";
  input Integer m1 "Order of the first harmonic";
  input Integer l2 "Degree of the second harmonic";
  input Integer m2 "Order of the second harmonic";
  input Integer l3 "Degree of the third harmonic";
  input Integer m3 "Order of the third harmonic";
  output Real G "Integral of Y_l1m1 Y_l2m2 Y_l3m3 dOmega (no complex conjugates)";
algorithm
  G := sqrt((2 * l1 + 1) * (2 * l2 + 1) * (2 * l3 + 1) / (4.0 * GeneralRelativity.Constants.pi))
    * wigner3j(l1, l2, l3, 0, 0, 0) * wigner3j(l1, l2, l3, m1, m2, m3);
  annotation(Documentation(info="<html>
<p><code>&int;Y<sub>l1m1</sub>Y<sub>l2m2</sub>Y<sub>l3m3</sub>d&Omega; =
sqrt((2l1+1)(2l2+1)(2l3+1)/(4&pi;)) (l1 l2 l3; 0 0 0)(l1 l2 l3; m1 m2 m3)</code>
with the complex harmonics of
<a href=\"modelica://GeneralRelativity.Harmonics.sphericalHarmonic\">sphericalHarmonic</a>.
It vanishes unless <code>m1 + m2 + m3 = 0</code>, the triangle condition holds and
<code>l1 + l2 + l3</code> is even (parity). For a conjugated harmonic use
<code>Y*<sub>lm</sub> = (-1)<sup>m</sup>Y<sub>l,-m</sub></code>, e.g. the coupling of
<code>Y<sub>l1m1</sub>Y<sub>l2m2</sub></code> onto <code>Y<sub>LM</sub></code> is
<code>(-1)<sup>M</sup> gaunt(l1, m1, l2, m2, L, -M)</code>.</p>
</html>"));
end gaunt;
