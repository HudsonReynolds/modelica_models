within;
package GeneralRelativity "Metric curvature, stress-energy, weak-field perturbation theory, geodesics and field solvers in the (-+++) convention"
  annotation(
    uses(Modelica(version="4.1.0"), LinearAlgebra),
    Documentation(info="<html>
<h4>Conventions</h4>
<p>Spacetime indices run <code>1..4</code> in Modelica arrays and correspond to
<code>&mu; = 0..3</code>: index 1 is time and 2..4 are space. The time
coordinate is <code>x[1] = c*t</code>, so every coordinate has the dimension of
the spatial chart (metres in Cartesian charts, metres and radians in
spherical charts). The signature is <b>(&minus;,+,+,+)</b>:
the Minkowski metric is <code>diag(-1, 1, 1, 1)</code> and a timelike
four-velocity satisfies <code>g(u,u) = -c&sup2;</code>.</p>
<p>Curvature follows Misner, Thorne and Wheeler and Carroll:</p>
<ul>
<li><code>&Gamma;<sup>i</sup><sub>jk</sub> = &frac12; g<sup>il</sup>(&part;<sub>j</sub>g<sub>lk</sub> + &part;<sub>k</sub>g<sub>lj</sub> &minus; &part;<sub>l</sub>g<sub>jk</sub>)</code>,
stored as <code>Gamma[i, j, k]</code>;</li>
<li><code>R<sup>&rho;</sup><sub>&sigma;&mu;&nu;</sub> = &part;<sub>&mu;</sub>&Gamma;<sup>&rho;</sup><sub>&nu;&sigma;</sub> &minus; &part;<sub>&nu;</sub>&Gamma;<sup>&rho;</sup><sub>&mu;&sigma;</sub> + &Gamma;<sup>&rho;</sup><sub>&mu;&lambda;</sub>&Gamma;<sup>&lambda;</sup><sub>&nu;&sigma;</sub> &minus; &Gamma;<sup>&rho;</sup><sub>&nu;&lambda;</sub>&Gamma;<sup>&lambda;</sup><sub>&mu;&sigma;</sub></code>,
so a sphere has positive curvature;</li>
<li><code>R<sub>&mu;&nu;</sub> = R<sup>&lambda;</sup><sub>&mu;&lambda;&nu;</sub></code>,
<code>R = g<sup>&mu;&nu;</sup>R<sub>&mu;&nu;</sub></code>,
<code>G<sub>&mu;&nu;</sub> = R<sub>&mu;&nu;</sub> &minus; &frac12; g<sub>&mu;&nu;</sub> R</code>;</li>
<li>field equations <code>G<sub>&mu;&nu;</sub> + &Lambda; g<sub>&mu;&nu;</sub> = (8&pi;G/c<sup>4</sup>) T<sub>&mu;&nu;</sub></code>.</li>
</ul>
<p>Partial derivatives are stored with the derivative index first:
<code>dg[k, i, j] = &part;<sub>k</sub> g<sub>ij</sub></code> and
<code>ddg[k, l, i, j] = &part;<sub>k</sub>&part;<sub>l</sub> g<sub>ij</sub></code>.
All tensors returned by this library carry lower indices unless the name says
otherwise (<code>riemannMixed</code> raises the first index).</p>

<h4>Metric functions</h4>
<p>A metric is a function with the signature of
<a href=\"modelica://GeneralRelativity.Interfaces.MetricFunction\">Interfaces.MetricFunction</a>:
coordinates and a parameter vector in, the covariant metric out. Curvature,
geodesic and visualization functions take the metric as a functional input, so
the same code evaluates
<a href=\"modelica://GeneralRelativity.Metrics.schwarzschild\">Schwarzschild</a>,
<a href=\"modelica://GeneralRelativity.Metrics.kerr\">Kerr</a>, the weak-field
<a href=\"modelica://GeneralRelativity.Metrics.earth\">Earth</a> and
<a href=\"modelica://GeneralRelativity.Metrics.cislunar\">cislunar</a> metrics or a
user function. Metric derivatives are central finite differences with a
caller-chosen step per coordinate
(<a href=\"modelica://GeneralRelativity.Metrics.curvatureAt\">Metrics.curvatureAt</a>).</p>

<h4>Choose a package</h4>
<table>
<tr><th>Task</th><th>Package</th></tr>
<tr><td>Christoffel, Riemann, Ricci, Einstein and stress-energy tensors from a metric and its derivatives</td>
<td><a href=\"modelica://GeneralRelativity.Tensors\">Tensors</a></td></tr>
<tr><td>Named metrics and numerical curvature at a point</td>
<td><a href=\"modelica://GeneralRelativity.Metrics\">Metrics</a></td></tr>
<tr><td>Vacuum, dust, perfect fluid and cosmological-constant sources</td>
<td><a href=\"modelica://GeneralRelativity.StressEnergy\">StressEnergy</a></td></tr>
<tr><td>Linearized gravity, Newtonian limit, J2 and two-body potentials, post-Newtonian motion</td>
<td><a href=\"modelica://GeneralRelativity.Perturbation\">Perturbation</a></td></tr>
<tr><td>Geodesic equation and a replaceable-metric geodesic model</td>
<td><a href=\"modelica://GeneralRelativity.Geodesics\">Geodesics</a></td></tr>
<tr><td>Finite-difference Poisson, spherical finite elements, method-of-lines wave propagation</td>
<td><a href=\"modelica://GeneralRelativity.Numerical\">Numerical</a></td></tr>
<tr><td>Embedding surfaces, curvature fields and the Three.js scene stream</td>
<td><a href=\"modelica://GeneralRelativity.Visualization\">Visualization</a></td></tr>
<tr><td>Schwarzschild orbit and curvature, Earth J2, cislunar, Poisson, star interior, wave pulse</td>
<td><a href=\"modelica://GeneralRelativity.Examples\">Examples</a></td></tr>
</table>

<h4>Units</h4>
<p>Every function takes <code>c</code> and <code>G</code> explicitly where they
enter. Pass SI values from
<a href=\"modelica://GeneralRelativity.Constants\">Constants</a> for Earth and
cislunar work, or <code>c = 1</code>, <code>G = 1</code> for geometric units,
where mass is a length (<code>GM/c&sup2;</code>) and the Schwarzschild radius is
<code>2M</code>. The weak-field examples use post-Newtonian equations of motion
rather than finite differences of a metric that differs from Minkowski by one
part in 10<sup>9</sup>; see
<a href=\"modelica://GeneralRelativity.Perturbation\">Perturbation</a> for why.</p>

<h4>Example</h4>
<pre>GeneralRelativity.Tensors.Curvature k =
  GeneralRelativity.Metrics.curvatureAt(
    GeneralRelativity.Metrics.schwarzschild, {0, 10, Modelica.Constants.pi/2, 0},
    {2.0}, {1e-3, 1e-3, 1e-3, 1e-3});
// k.ricci is zero, k.kretschmann is 12*rs^2/r^6 = 48e-6</pre>
<p><a href=\"modelica://GeneralRelativity.Examples.SchwarzschildOrbit\">Examples.SchwarzschildOrbit</a>
integrates a precessing orbit and checks the conserved energy and angular momentum;
<a href=\"modelica://GeneralRelativity.Examples.SchwarzschildOrbitScene\">SchwarzschildOrbitScene</a>
adds the scene stream for <code>Visualization/threejs/index.html</code>.</p>
</html>"));
end GeneralRelativity;
