within GeneralRelativity;
package Metrics "Named metrics, finite-difference metric derivatives and curvature at a point"
  annotation(Documentation(info="<html>
<h4>Charts and parameter vectors</h4>
<table>
<tr><th>Metric</th><th>Coordinates</th><th><code>p</code></th></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Metrics.minkowski\">minkowski</a></td><td>{ct, x, y, z}</td><td>any</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Metrics.schwarzschild\">schwarzschild</a></td><td>{ct, r, &theta;, &phi;}</td><td>{r<sub>s</sub>}</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Metrics.kerr\">kerr</a></td><td>{ct, r, &theta;, &phi;} Boyer-Lindquist</td><td>{r<sub>s</sub>, a}</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Metrics.weakFieldPointMass\">weakFieldPointMass</a></td><td>{ct, x, y, z}</td><td>{GM, c}</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Metrics.earth\">earth</a></td><td>{ct, x, y, z}, z along the spin axis</td><td>{GM, R, J2, c, J}</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Metrics.cislunar\">cislunar</a></td><td>{ct, x, y, z} Earth-centred inertial</td><td>{GM<sub>E</sub>, GM<sub>M</sub>, d, &Omega;, c, &phi;<sub>0</sub>}</td></tr>
</table>
<h4>Finite differences</h4>
<p><a href=\"modelica://GeneralRelativity.Metrics.derivative\">derivative</a> and
<a href=\"modelica://GeneralRelativity.Metrics.secondDerivative\">secondDerivative</a>
use central differences with a step per coordinate. Truncation error is
<code>O(h&sup2;)</code> and round-off is <code>O(&epsilon;/h&sup2;)</code> for the second
derivative, so a step near <code>10<sup>-3</sup></code> to <code>10<sup>-4</sup></code>
of the coordinate scale gives six to eight significant figures for
<code>O(1)</code> metrics. For a metric within <code>10<sup>-9</sup></code> of
Minkowski the curvature signal is below that floor; use the linearized tools in
<a href=\"modelica://GeneralRelativity.Perturbation\">Perturbation</a> with
analytic potentials instead.</p>
<p>Angular coordinates need radian-sized steps, and the step must stay away from
chart singularities: the horizon <code>r = r<sub>s</sub></code> and the axis
<code>&theta; = 0, &pi;</code>.</p>
</html>"));
end Metrics;
