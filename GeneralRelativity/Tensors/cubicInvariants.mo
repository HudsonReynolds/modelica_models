within GeneralRelativity.Tensors;
function cubicInvariants "Basis of six cubic curvature invariants {R^3, R Ric^2, Ric^3, R_ab R_cd R^acbd, R K, Riem^3}, in any dimension"
  input Real gInv[:, size(gInv, 1)] "Contravariant metric";
  input Real RLower[size(gInv, 1), size(gInv, 1), size(gInv, 1), size(gInv, 1)] "R_abcd";
  input Real Ric[size(gInv, 1), size(gInv, 1)] "Ricci tensor R_ab";
  input Real R "Ricci scalar";
  output Real I[6] "Cubic invariants in the basis documented below";
protected
  Integer n = size(gInv, 1);
  Real mixed[size(gInv, 1), size(gInv, 1)] "R^a_b";
  Real RicUp[size(gInv, 1), size(gInv, 1)] "R^ab";
  Real M[size(gInv, 1) * size(gInv, 1), size(gInv, 1) * size(gInv, 1)] "R_ab^cd";
  Real M2[size(gInv, 1) * size(gInv, 1), size(gInv, 1) * size(gInv, 1)];
  Real Ric2;
  Real K;
algorithm
  mixed := gInv * Ric;
  RicUp := mixed * gInv;
  M := bivectorMatrix(gInv, RLower);
  M2 := M * M;
  Ric2 := sum(mixed[a, b] * mixed[b, a] for a in 1:n, b in 1:n);
  K := sum(M2[i, i] for i in 1:n * n);
  I[1] := R^3;
  I[2] := R * Ric2;
  I[3] := sum(mixed[a, b] * mixed[b, c] * mixed[c, a] for a in 1:n, b in 1:n, c in 1:n);
  I[4] := sum(RicUp[a, b] * RicUp[c, d] * RLower[a, c, b, d] for a in 1:n, b in 1:n, c in 1:n, d in 1:n);
  I[5] := R * K;
  I[6] := sum(M2[i, j] * M[j, i] for i in 1:n * n, j in 1:n * n);
  annotation(Documentation(info="<html>
<p>The six algebraically independent cubic scalars that can be formed from the Riemann tensor
without derivatives (generic dimension; the second Riemann-cubed contraction
<code>R<sub>a&nbsp;b</sub><sup>&nbsp;c&nbsp;d</sup>R<sub>c&nbsp;d</sub><sup>&nbsp;e&nbsp;f</sup>R<sub>e&nbsp;f</sub><sup>&nbsp;a&nbsp;b</sup></code> with interleaved
pairs is left out):</p>
<table>
<tr><th>Index</th><th>Invariant</th></tr>
<tr><td>1</td><td><code>R&sup3;</code></td></tr>
<tr><td>2</td><td><code>R R<sub>ab</sub>R<sup>ab</sup></code></td></tr>
<tr><td>3</td><td><code>R<sub>a</sub><sup>b</sup>R<sub>b</sub><sup>c</sup>R<sub>c</sub><sup>a</sup></code></td></tr>
<tr><td>4</td><td><code>R<sub>ab</sub>R<sub>cd</sub>R<sup>acbd</sup></code></td></tr>
<tr><td>5</td><td><code>R R<sub>abcd</sub>R<sup>abcd</sup></code></td></tr>
<tr><td>6</td><td><code>R<sub>ab</sub><sup>cd</sup>R<sub>cd</sub><sup>ef</sup>R<sub>ef</sub><sup>ab</sup></code></td></tr>
</table>
<p>For an Einstein space <code>R<sub>ab</sub> = &Lambda;g<sub>ab</sub></code> in four dimensions
the first four are <code>64&Lambda;&sup3;, 16&Lambda;&sup3;, 4&Lambda;&sup3;, 4&Lambda;&sup3;</code>.
In vacuum only entry 6 survives; for Schwarzschild it is <code>12 r<sub>s</sub>&sup3;/r<sup>9</sup></code>.
The Weyl analogue is <a href=\"modelica://GeneralRelativity.Tensors.weylCubed\">weylCubed</a>.</p>
</html>"));
end cubicInvariants;
