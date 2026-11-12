within GeneralRelativity;
package Interfaces "Function signatures passed as functional inputs"
  annotation(Documentation(info="<html>
<p>Curvature, geodesic and visualization routines take a metric as a
<i>functional input</i>, Modelica specification section 12.4.2. Any function
whose inputs and outputs match
<a href=\"modelica://GeneralRelativity.Interfaces.MetricFunction\">MetricFunction</a>
can be passed by name, for example
<code>Metrics.curvatureAt(GeneralRelativity.Metrics.kerr, x, {rs, a}, step)</code>.
Metric parameters travel in a plain vector so a model can bind them to its own
parameters; each metric documents its vector layout.</p>
</html>"));
end Interfaces;
