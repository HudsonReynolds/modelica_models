within GeneralRelativity.Frames;
function weylInvariants "Complex Weyl invariants I = Psi0 Psi4 - 4 Psi1 Psi3 + 3 Psi2^2 and J = det of the Psi Hankel matrix"
  input Real psiRe[5] "Real parts of Psi0..Psi4";
  input Real psiIm[5] "Imaginary parts of Psi0..Psi4";
  output Real IRe "Re I";
  output Real IIm "Im I";
  output Real JRe "Re J";
  output Real JIm "Im J";
protected
  // A complex number z = x + i y is the real matrix [x, -y; y, x]; products are matrix products.
  Real P[5, 2, 2];
  Real Ic[2, 2];
  Real Jc[2, 2];
algorithm
  for k in 1:5 loop
    P[k, :, :] := [psiRe[k], -psiIm[k]; psiIm[k], psiRe[k]];
  end for;
  Ic := P[1, :, :] * P[5, :, :] - 4.0 * P[2, :, :] * P[4, :, :] + 3.0 * P[3, :, :] * P[3, :, :];
  // det [Psi4, Psi3, Psi2; Psi3, Psi2, Psi1; Psi2, Psi1, Psi0] by the first row.
  Jc := P[5, :, :] * (P[3, :, :] * P[1, :, :] - P[2, :, :] * P[2, :, :])
    - P[4, :, :] * (P[4, :, :] * P[1, :, :] - P[2, :, :] * P[3, :, :])
    + P[3, :, :] * (P[4, :, :] * P[2, :, :] - P[3, :, :] * P[3, :, :]);
  IRe := Ic[1, 1];
  IIm := Ic[2, 1];
  JRe := Jc[1, 1];
  JIm := Jc[2, 1];
  annotation(Documentation(info="<html>
<p><code>I = &Psi;<sub>0</sub>&Psi;<sub>4</sub> &minus; 4&Psi;<sub>1</sub>&Psi;<sub>3</sub> + 3&Psi;<sub>2</sub>&sup2;</code> and
<code>J = det[[&Psi;<sub>4</sub>, &Psi;<sub>3</sub>, &Psi;<sub>2</sub>], [&Psi;<sub>3</sub>, &Psi;<sub>2</sub>, &Psi;<sub>1</sub>], [&Psi;<sub>2</sub>, &Psi;<sub>1</sub>, &Psi;<sub>0</sub>]]</code>
are unchanged by null rotations and boosts of the tetrad. The Weyl tensor is algebraically special
(types II, D, III, N, O) exactly when <code>I&sup3; = 27J&sup2;</code>. With the definitions of
<a href=\"modelica://GeneralRelativity.Frames.weylScalars\">weylScalars</a>, <code>16 Re I = C<sub>abcd</sub>C<sup>abcd</sup></code>.
For a type D field with only <code>&Psi;<sub>2</sub></code> nonzero, <code>I = 3&Psi;<sub>2</sub>&sup2;</code> and
<code>J = &minus;&Psi;<sub>2</sub>&sup3;</code>.</p>
</html>"));
end weylInvariants;
