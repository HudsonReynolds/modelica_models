within GeneralRelativity.Tensors;
function schouten "Schouten tensor P_ab = (R_ab - R g_ab / (2 (n - 1))) / (n - 2), in any dimension n >= 3"
  input Real g[:, size(g, 1)] "Covariant metric";
  input Real Ric[size(g, 1), size(g, 1)] "Ricci tensor R_ab";
  input Real R "Ricci scalar";
  output Real P[size(g, 1), size(g, 1)] "P_ab";
protected
  Integer n = size(g, 1);
algorithm
  assert(n >= 3, "schouten: the Schouten tensor needs dimension 3 or more");
  P := (Ric - R / (2.0 * (n - 1)) * g) / (n - 2);
  annotation(Documentation(info="<html>
<p>The Riemann tensor splits as <code>R<sub>abcd</sub> = C<sub>abcd</sub> + (g &#8743; P)<sub>abcd</sub></code>
with <code>(g &#8743; P)<sub>abcd</sub> = g<sub>ac</sub>P<sub>bd</sub> &minus; g<sub>ad</sub>P<sub>bc</sub>
&minus; g<sub>bc</sub>P<sub>ad</sub> + g<sub>bd</sub>P<sub>ac</sub></code>; its trace is
<code>P = R/(2(n&minus;1))</code>. Under a conformal rescaling it transforms by second derivatives of
the conformal factor only, which is why the
<a href=\"modelica://GeneralRelativity.Metrics.cottonAt\">Cotton tensor</a> built from it measures
conformal flatness in three dimensions and the source of the Weyl tensor in four.</p>
</html>"));
end schouten;
