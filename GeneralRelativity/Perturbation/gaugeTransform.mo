within GeneralRelativity.Perturbation;
function gaugeTransform "Linearized gauge transformation h_ab -> h_ab - d_a xi_b - d_b xi_a on a flat background"
  input Real h[4, 4] "Perturbation h_ab (or any derivative d_k d_l h_ab of it)";
  input Real dxi[4, 4] "dxi[a, b] = d_a xi_b (or the matching derivative d_k d_l d_a xi_b)";
  output Real hNew[4, 4] "Transformed perturbation h_ab - d_a xi_b - d_b xi_a";
algorithm
  hNew := h - dxi - transpose(dxi);
  annotation(Documentation(info="<html>
<p>An infinitesimal coordinate change <code>x<sup>a</sup> &rarr; x<sup>a</sup> + &xi;<sup>a</sup></code>
changes the perturbation of <code>g = &eta; + h</code> by the Lie derivative of
<code>&eta;</code>:
<code>h<sub>ab</sub> &rarr; h<sub>ab</sub> - &part;<sub>a</sub>&xi;<sub>b</sub> - &part;<sub>b</sub>&xi;<sub>a</sub></code>
with <code>&xi;<sub>b</sub> = &eta;<sub>bc</sub>&xi;<sup>c</sup></code>. The linearized
Riemann tensor
(<a href=\"modelica://GeneralRelativity.Perturbation.linearizedRiemann\">linearizedRiemann</a>)
is invariant, so curvature, tidal fields and gravitational-wave strain measured by
<code>R<sub>0i0j</sub></code> are gauge independent, while <code>h</code> itself is not.
A Killing vector of Minkowski space (translation, rotation or boost) has
antisymmetric <code>&part;<sub>a</sub>&xi;<sub>b</sub></code> and leaves <code>h</code> unchanged.</p>
<p>Because partial derivatives commute, the same function transforms derivatives of
<code>h</code>: pass <code>ddh[k, l, :, :]</code> and
<code>&part;<sub>k</sub>&part;<sub>l</sub>&part;<sub>a</sub>&xi;<sub>b</sub></code> to obtain the
transformed <code>ddh[k, l, :, :]</code>.</p>
</html>"));
end gaugeTransform;
