within GeneralRelativity.Foliation;
function admAt "3+1 split, extrinsic curvature, 3-curvature and Gauss-Codazzi residuals of a metric function at one event"
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric function";
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Metric parameters";
  input Real step[4] "Finite-difference step per coordinate";
  output ADM s "3+1 quantities and constraint residuals";
protected
  Real g[4, 4];
  Real dg[4, 4, 4];
  Real ddg[4, 4, 4, 4];
algorithm
  g := metric(x, p);
  dg := GeneralRelativity.Metrics.derivative(metric, x, p, step);
  ddg := GeneralRelativity.Metrics.secondDerivative(metric, x, p, step);
  s := admFromJet(x, g, dg, ddg);
  annotation(Documentation(info="<html>
<p>Central differences as in
<a href=\"modelica://GeneralRelativity.Metrics.curvatureAt\">Metrics.curvatureAt</a>
(41 metric evaluations), then
<a href=\"modelica://GeneralRelativity.Foliation.admFromJet\">admFromJet</a>.
The time derivative of <code>&gamma;<sub>ij</sub></code> uses <code>step[1]</code>
in <code>x<sup>0</sup> = ct</code>. The spatial derivatives of <code>K<sub>ij</sub></code>
needed by the momentum constraint come from the chain rule on the jet, not
from differencing <code>K</code> at neighbouring events.</p>
<pre>GeneralRelativity.Foliation.ADM s = GeneralRelativity.Foliation.admAt(
  GeneralRelativity.Metrics.schwarzschild, {0, 10, 1.2, 0.3}, {2}, fill(1e-3, 4));
// time-symmetric slice: s.K = 0, s.ricciScalar3 = 0, s.alpha = sqrt(1 - 2/10)</pre>
</html>"));
end admAt;
