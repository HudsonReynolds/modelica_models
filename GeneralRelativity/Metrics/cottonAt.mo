within GeneralRelativity.Metrics;
function cottonAt "Cotton tensor C_abc = nabla_c P_ab - nabla_b P_ac of a metric function at one event, by finite differences of the Schouten tensor"
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric function";
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Metric parameters";
  input Real step[4] "Finite-difference step per coordinate (used for the curvature and for its derivative)";
  output Real C[4, 4, 4] "C[a, b, c] = C_abc";
protected
  GeneralRelativity.Tensors.Curvature k0;
  GeneralRelativity.Tensors.Curvature kPlus;
  GeneralRelativity.Tensors.Curvature kMinus;
  Real P[4, 4] "Schouten tensor at x";
  Real dP[4, 4, 4] "dP[c, a, b] = d_c P_ab";
  Real nablaP[4, 4, 4] "nablaP[c, a, b] = nabla_c P_ab";
  Real xs[4];
algorithm
  k0 := curvatureAt(metric, x, p, step);
  P := GeneralRelativity.Tensors.schouten(k0.g, k0.ricci, k0.ricciScalar);
  for c in 1:4 loop
    xs := x;
    xs[c] := x[c] + step[c];
    kPlus := curvatureAt(metric, xs, p, step);
    xs[c] := x[c] - step[c];
    kMinus := curvatureAt(metric, xs, p, step);
    dP[c, :, :] := (GeneralRelativity.Tensors.schouten(kPlus.g, kPlus.ricci, kPlus.ricciScalar)
      - GeneralRelativity.Tensors.schouten(kMinus.g, kMinus.ricci, kMinus.ricciScalar)) / (2.0 * step[c]);
  end for;
  for c in 1:4 loop
    for a in 1:4 loop
      for b in 1:4 loop
        nablaP[c, a, b] := dP[c, a, b] - sum(k0.Gamma[e, c, a] * P[e, b] + k0.Gamma[e, c, b] * P[a, e] for e in 1:4);
      end for;
    end for;
  end for;
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        C[a, b, c] := nablaP[c, a, b] - nablaP[b, a, c];
      end for;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>With the <a href=\"modelica://GeneralRelativity.Tensors.schouten\">Schouten tensor</a>
<code>P<sub>ab</sub></code>, <code>C<sub>abc</sub> = &nabla;<sub>c</sub>P<sub>ab</sub> &minus; &nabla;<sub>b</sub>P<sub>ac</sub></code>.
It is antisymmetric in <code>b, c</code>, its cyclic sum vanishes and, by the contracted Bianchi
identity, both traces <code>g<sup>ab</sup>C<sub>abc</sub></code> and <code>g<sup>ac</sup>C<sub>abc</sub></code> vanish.
In four dimensions it is the divergence of the Weyl tensor,
<code>&nabla;<sup>d</sup>C<sub>dcab</sub> = C<sub>abc</sub></code> up to sign convention, so it vanishes for
conformally flat metrics and for Einstein spaces (where <code>P &prop; g</code>). In three dimensions it
vanishes exactly when the metric is conformally flat.</p>
<p>Cost: nine calls of <a href=\"modelica://GeneralRelativity.Metrics.curvatureAt\">curvatureAt</a>
(the event and two neighbours per coordinate), i.e. 297 metric evaluations. The nested central
differences are second-order accurate in <code>step</code>; with a step of <code>10<sup>&minus;3</sup></code>
of the coordinate scale the round-off is about <code>10<sup>&minus;7</sup></code> in units of the metric over
the coordinate scale cubed.</p>
</html>"));
end cottonAt;
