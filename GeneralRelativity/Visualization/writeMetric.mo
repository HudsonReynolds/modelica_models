within GeneralRelativity.Visualization;
impure function writeMetric "Declare the scene's metric by family name and parameters, for viewer styles that evaluate it anywhere (the ray tracer)"
  input String fileName;
  input String name "Metric family known to the viewer: minkowski, schwarzschild, kerr, reissnerNordstrom or kerrNewman";
  input Real p[:] "Parameters in the layout of the Modelica function of the same name in Metrics";
  input String chart = "spherical" "Chart of the scene samples";
algorithm
  Modelica.Utilities.Streams.print("{\"type\":\"metric\",\"name\":\"" + name + "\",\"chart\":\"" + chart
    + "\",\"p\":" + jsonVector(p) + "}", fileName);
  annotation(Documentation(info="<html>
<p>The viewer re-implements these families in GLSL in Kerr-Schild form,
<code>g = &eta; + f k&otimes;k</code>, which is regular across the horizon; the same formula is
<a href=\"modelica://GeneralRelativity.Metrics.kerrSchild\">Metrics.kerrSchild</a> and is tested against
the Boyer-Lindquist metrics. Scenes with other metrics simply omit this record; the slice-based styles
still work.</p>
</html>"));
end writeMetric;
