within GeneralRelativity.Symmetry;
function lieDerivativeMetric "Lie derivative of the metric along a vector field, (L_xi g)_ab"
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric function";
  input GeneralRelativity.Interfaces.VectorFieldFunction field "Vector field xi";
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Metric parameters";
  input Real q[:] "Field parameters";
  input Real step[4] "Finite-difference step per coordinate";
  output Real L[4, 4] "(L_xi g)_ab = xi^c d_c g_ab + g_cb d_a xi^c + g_ac d_b xi^c";
  output Real scale "Largest sum of the magnitudes of the terms in one component; the natural size of L";
protected
  Real g[4, 4];
  Real dg[4, 4, 4];
  Real xi[4];
  Real dxi[4, 4];
  Real magnitude;
algorithm
  g := metric(x, p);
  dg := GeneralRelativity.Metrics.derivative(metric, x, p, step);
  xi := field(x, q);
  dxi := fieldDerivative(field, x, q, step);
  scale := 0.0;
  for a in 1:4 loop
    for b in 1:4 loop
      L[a, b] := 0.0;
      magnitude := 0.0;
      for c in 1:4 loop
        L[a, b] := L[a, b] + xi[c] * dg[c, a, b] + g[c, b] * dxi[a, c] + g[a, c] * dxi[b, c];
        magnitude := magnitude + abs(xi[c] * dg[c, a, b]) + abs(g[c, b] * dxi[a, c]) + abs(g[a, c] * dxi[b, c]);
      end for;
      scale := max(scale, magnitude);
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>Equivalent to <code>&nabla;<sub>a</sub>&xi;<sub>b</sub> + &nabla;<sub>b</sub>&xi;<sub>a</sub></code>
but needs no Christoffel symbols. Metric and field derivatives are central
differences (<a href=\"modelica://GeneralRelativity.Metrics.derivative\">Metrics.derivative</a>,
<a href=\"modelica://GeneralRelativity.Symmetry.fieldDerivative\">fieldDerivative</a>).
<code>scale</code> lets callers judge cancellation; see
<a href=\"modelica://GeneralRelativity.Symmetry.killingResidual\">killingResidual</a>.</p>
</html>"));
end lieDerivativeMetric;
