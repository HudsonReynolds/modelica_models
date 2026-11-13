within GeneralRelativity.Metrics;
function flrwCurved "FLRW metric with spatial curvature k in comoving spherical coordinates {ct, r, theta, phi}; p = {ct0, n, k}"
  extends GeneralRelativity.Interfaces.MetricFunction;
protected
  Real ct0 = p[1] "Reference time c*t0 at which a = 1";
  Real n = p[2] "Power-law exponent of a = (ct/ct0)^n";
  Real k = p[3] "Spatial curvature (1/length^2): positive closed, negative open";
  Real r = x[2];
  Real a = (x[1] / ct0)^n;
algorithm
  g := zeros(4, 4);
  g[1, 1] := -1.0;
  g[2, 2] := a^2 / (1.0 - k * r^2);
  g[3, 3] := a^2 * r^2;
  g[4, 4] := a^2 * r^2 * sin(x[3])^2;
  annotation(Documentation(info="<html>
<p><code>ds&sup2; = -c&sup2;dt&sup2; + a&sup2;[dr&sup2;/(1 - k r&sup2;) + r&sup2;(d&theta;&sup2; + sin&sup2;&theta; d&phi;&sup2;)]</code>
with the power-law scale factor <code>a = (ct/ct0)<sup>n</sup></code>. The comoving radius
<code>r</code> has the dimension of length, <code>k</code> of 1/length&sup2;, and the
chart needs <code>k r&sup2; &lt; 1</code>. The spatial slices have Ricci scalar
<code>6k/a&sup2;</code>; the Friedmann constraint is
<code>G<sub>00</sub> = 3(H&sup2; + k/a&sup2;)</code> with <code>H = (da/d(ct))/a</code>,
and the spatial equations give
<code>G<sup>r</sup><sub>r</sub> = -(2&auml;/a + H&sup2; + k/a&sup2;)</code>.
A power law solves them with a single perfect fluid only when <code>k = 0</code>
(<a href=\"modelica://GeneralRelativity.Metrics.flrw\">flrw</a>); otherwise read the
implied source off the Einstein tensor.</p>
</html>"));
end flrwCurved;
