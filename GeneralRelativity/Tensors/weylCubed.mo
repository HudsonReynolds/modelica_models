within GeneralRelativity.Tensors;
function weylCubed "Cubic invariant C_ab^cd C_cd^ef C_ef^ab (Goroff-Sagnotti); also valid for any tensor with Riemann symmetries"
  input Real gInv[:, size(gInv, 1)] "Contravariant metric";
  input Real C[size(gInv, 1), size(gInv, 1), size(gInv, 1), size(gInv, 1)] "C_abcd (Weyl, or Riemann for R_ab^cd R_cd^ef R_ef^ab)";
  output Real C3 "C_ab^cd C_cd^ef C_ef^ab";
protected
  Integer n = size(gInv, 1);
  Real M[size(gInv, 1) * size(gInv, 1), size(gInv, 1) * size(gInv, 1)] "M[(ab), (cd)] = C_ab^cd";
  Real M2[size(gInv, 1) * size(gInv, 1), size(gInv, 1) * size(gInv, 1)];
algorithm
  M := bivectorMatrix(gInv, C);
  M2 := M * M;
  C3 := sum(M2[i, j] * M[j, i] for i in 1:n * n, j in 1:n * n);
  annotation(Documentation(info="<html>
<p>The tensor is read as a linear map <code>M</code> on pairs of indices,
<code>M<sub>(ab)(cd)</sub> = C<sub>ab</sub><sup>cd</sup></code>, and the invariant is
<code>tr M&sup3;</code>. This is the two-loop counterterm of pure Einstein gravity found by
Goroff and Sagnotti. For Schwarzschild, where the Weyl tensor equals the Riemann tensor,
the value is <code>96 M&sup3;/r<sup>9</sup> = 12 r<sub>s</sub>&sup3;/r<sup>9</sup></code>.</p>
</html>"));
end weylCubed;
