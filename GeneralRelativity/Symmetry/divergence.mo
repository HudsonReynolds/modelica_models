within GeneralRelativity.Symmetry;
function divergence "Covariant divergence nabla_a xi^a of a vector field"
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric function";
  input GeneralRelativity.Interfaces.VectorFieldFunction field "Vector field xi";
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Metric parameters";
  input Real q[:] "Field parameters";
  input Real step[4] "Finite-difference step per coordinate";
  output Real div "nabla_a xi^a = d_a xi^a + (1/2) g^bc (xi^a d_a g_bc)";
protected
  Real g[4, 4];
  Real gInv[4, 4];
  Boolean accepted;
  Real dg[4, 4, 4];
  Real xi[4];
  Real dxi[4, 4];
algorithm
  g := metric(x, p);
  (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(g);
  dg := GeneralRelativity.Metrics.derivative(metric, x, p, step);
  xi := field(x, q);
  dxi := fieldDerivative(field, x, q, step);
  div := sum(dxi[a, a] for a in 1:4);
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        div := div + 0.5 * xi[a] * gInv[b, c] * dg[a, b, c];
      end for;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>Uses <code>&Gamma;<sup>a</sup><sub>ab</sub> = &part;<sub>b</sub> ln&radic;|g| = &frac12;g<sup>ac</sup>&part;<sub>b</sub>g<sub>ac</sub></code>.
For <code>&part;<sub>r</sub></code> in a spherical chart of Schwarzschild the
result is <code>2/r</code>; for the dilation in Minkowski it is 4.</p>
</html>"));
end divergence;
