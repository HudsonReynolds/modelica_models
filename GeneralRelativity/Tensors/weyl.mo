within GeneralRelativity.Tensors;
function weyl "Weyl conformal tensor C_abcd: the trace-free part of the Riemann tensor, in any dimension n >= 3"
  input Real g[:, size(g, 1)] "Covariant metric";
  input Real RLower[size(g, 1), size(g, 1), size(g, 1), size(g, 1)] "R_abcd";
  input Real Ric[size(g, 1), size(g, 1)] "Ricci tensor R_ab";
  input Real R "Ricci scalar";
  output Real C[size(g, 1), size(g, 1), size(g, 1), size(g, 1)] "C_abcd";
protected
  Integer n = size(g, 1);
  Real P[size(g, 1), size(g, 1)] "Schouten tensor";
algorithm
  assert(n >= 3, "weyl: the Weyl tensor needs dimension 3 or more");
  P := schouten(g, Ric, R);
  for a in 1:n loop
    for b in 1:n loop
      for c in 1:n loop
        for d in 1:n loop
          C[a, b, c, d] := RLower[a, b, c, d]
            - (g[a, c] * P[b, d] - g[a, d] * P[b, c] - g[b, c] * P[a, d] + g[b, d] * P[a, c]);
        end for;
      end for;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p><code>C<sub>abcd</sub> = R<sub>abcd</sub> &minus; (g<sub>ac</sub>P<sub>bd</sub> &minus; g<sub>ad</sub>P<sub>bc</sub>
&minus; g<sub>bc</sub>P<sub>ad</sub> + g<sub>bd</sub>P<sub>ac</sub>)</code> with the
<a href=\"modelica://GeneralRelativity.Tensors.schouten\">Schouten tensor</a> <code>P</code>.
Expanded, this is the familiar
<code>R<sub>abcd</sub> &minus; (g<sub>ac</sub>R<sub>bd</sub> &minus; g<sub>ad</sub>R<sub>bc</sub> &minus; g<sub>bc</sub>R<sub>ad</sub> + g<sub>bd</sub>R<sub>ac</sub>)/(n&minus;2)
+ R (g<sub>ac</sub>g<sub>bd</sub> &minus; g<sub>ad</sub>g<sub>bc</sub>)/((n&minus;1)(n&minus;2))</code>.
<code>C</code> has the Riemann symmetries and every trace <code>g<sup>ac</sup>C<sub>abcd</sub></code> vanishes.
In vacuum it equals the Riemann tensor; it vanishes identically for <code>n = 3</code> and for
conformally flat metrics <code>g = &Omega;&sup2;&eta;</code> when <code>n &ge; 4</code>.</p>
</html>"));
end weyl;
