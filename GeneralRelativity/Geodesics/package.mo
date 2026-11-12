within GeneralRelativity;
package Geodesics "Geodesic equation, four-velocity normalization, Killing constants and a replaceable-metric geodesic model"
  annotation(Documentation(info="<html>
<p>The geodesic equation in an affine parameter <code>&lambda;</code> is
<code>d&sup2;x<sup>i</sup>/d&lambda;&sup2; = -&Gamma;<sup>i</sup><sub>jk</sub> u<sup>j</sup>u<sup>k</sup></code>,
<code>u = dx/d&lambda;</code>. For massive particles <code>&lambda;</code> is proper time and
<code>g(u,u) = -c&sup2;</code>; for light <code>g(u,u) = 0</code>.
<a href=\"modelica://GeneralRelativity.Geodesics.Geodesic\">Geodesic</a> integrates it
with Modelica <code>time</code> as the affine parameter and any
<a href=\"modelica://GeneralRelativity.Interfaces.MetricFunction\">MetricFunction</a>
as the replaceable metric.</p>
<p>Stationary and axisymmetric metrics have the Killing vectors
<code>&part;<sub>t</sub></code> and <code>&part;<sub>&phi;</sub></code>;
<a href=\"modelica://GeneralRelativity.Geodesics.killingConserved\">killingConserved</a>
returns the corresponding constants <code>g<sub>ab</sub>&xi;<sup>a</sup>u<sup>b</sup></code>,
which the examples use as integration checks.</p>
</html>"));
end Geodesics;
