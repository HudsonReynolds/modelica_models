within GeneralRelativity.Geodesics;
function geodesicAcceleration "Geodesic acceleration of a metric function at x, with the metric and Christoffel symbols used"
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric function";
  input Real x[4] "Coordinates";
  input Real p[:] "Metric parameters";
  input Real step[4] "Finite-difference step per coordinate";
  input Real u[4] "Four-velocity dx/d lambda";
  output Real a[4] "d^2 x / d lambda^2";
  output Real g[4, 4] "Covariant metric at x";
  output Real Gamma[4, 4, 4] "Christoffel symbols at x";
protected
  Real gInv[4, 4];
  Boolean accepted;
  Real dg[4, 4, 4];
algorithm
  g := metric(x, p);
  (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(g);
  dg := GeneralRelativity.Metrics.derivative(metric, x, p, step);
  Gamma := GeneralRelativity.Tensors.christoffel(gInv, dg);
  a := acceleration(Gamma, u);
end geodesicAcceleration;
