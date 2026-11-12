within GeneralRelativity.Metrics;
function curvatureAt "Metric, Christoffel symbols and all curvature tensors of a metric function at one event"
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric function";
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Metric parameters";
  input Real step[4] "Finite-difference step per coordinate";
  output GeneralRelativity.Tensors.Curvature k;
protected
  Real dg[4, 4, 4];
  Real ddg[4, 4, 4, 4];
  Boolean accepted;
algorithm
  k.x := x;
  k.g := metric(x, p);
  (k.gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(k.g);
  k.valid := if accepted then 1.0 else 0.0;
  dg := derivative(metric, x, p, step);
  ddg := secondDerivative(metric, x, p, step);
  k.Gamma := GeneralRelativity.Tensors.christoffel(k.gInv, dg);
  k.riemann := GeneralRelativity.Tensors.riemannLower(k.g, ddg, k.Gamma);
  k.ricci := GeneralRelativity.Tensors.ricci(k.gInv, k.riemann);
  k.ricciScalar := GeneralRelativity.Tensors.ricciScalar(k.gInv, k.ricci);
  k.einstein := GeneralRelativity.Tensors.einstein(k.g, k.ricci, k.ricciScalar);
  k.kretschmann := GeneralRelativity.Tensors.kretschmann(k.gInv, k.riemann);
  annotation(Documentation(info="<html>
<p>Thirty-three metric evaluations: one at <code>x</code>, eight for the diagonal
second differences (which also give the first differences) and twenty-four for
mixed second differences. The stress-energy implied by the metric is
<code>Tensors.stressEnergyFromEinstein(k.einstein, k.g, G, c)</code>.</p>
</html>"));
end curvatureAt;
