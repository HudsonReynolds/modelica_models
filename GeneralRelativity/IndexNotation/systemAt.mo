within GeneralRelativity.IndexNotation;
function systemAt "System of curvature tensors of a metric function at one event, by central differences"
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric function";
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Metric parameters";
  input Real step[4] "Finite-difference step per coordinate";
  input Real gravitationalConstant(unit = "m3/(kg.s2)") = 1.0 "Newton's constant for T; 1 in geometric units";
  input Real c(unit = "m/s") = 1.0 "Speed of light for T; 1 in geometric units";
  input Real Lambda(unit = "1/m2") = 0.0 "Cosmological constant for T";
  output TensorSystem sys;
protected
  GeneralRelativity.Tensors.Curvature k;
algorithm
  k := GeneralRelativity.Metrics.curvatureAt(metric, x, p, step);
  assert(k.valid > 0.5, "IndexNotation.systemAt: the metric is singular at the requested event");
  sys := systemFromCurvature(k, gravitationalConstant, c, Lambda);
  annotation(Documentation(info="<html>
<p>Calls <a href=\"modelica://GeneralRelativity.Metrics.curvatureAt\">Metrics.curvatureAt</a>
and then <a href=\"modelica://GeneralRelativity.IndexNotation.systemFromCurvature\">systemFromCurvature</a>.</p>
<pre>sys = GeneralRelativity.IndexNotation.systemAt(
  GeneralRelativity.Metrics.kerr, {0, 6, pi/2, 0}, {2.0, 0.6}, fill(1e-3, 4));
GeneralRelativity.IndexNotation.scalar(sys, \"C_{abcd} C^{abcd}\");</pre>
</html>"));
end systemAt;
