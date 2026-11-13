within GeneralRelativity.Connections;
function leviCivita "Levi-Civita connection of a metric function at x, by central differences"
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric function";
  input Real x[4] "Coordinates";
  input Real p[:] "Metric parameters";
  input Real step[4] "Finite-difference step per coordinate";
  output Real Gamma[4, 4, 4] "Christoffel symbols Gamma^i_jk";
protected
  Real g[4, 4];
  Real gInv[4, 4];
  Boolean accepted;
  Real dg[4, 4, 4];
algorithm
  g := metric(x, p);
  (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(g);
  dg := GeneralRelativity.Metrics.derivative(metric, x, p, step);
  Gamma := GeneralRelativity.Tensors.christoffel(gInv, dg);
end leviCivita;
