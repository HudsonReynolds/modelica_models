within GeneralRelativity;
package Perturbation "Linearized gravity, the Newtonian and post-Newtonian limits, and weak-field potentials"
  annotation(Documentation(info="<html>
<h4>Linearized gravity</h4>
<p>Write <code>g = &eta; + h</code> with <code>|h| &lt;&lt; 1</code>. The linear Riemann
tensor is
<code>R<sub>abcd</sub> = &frac12;(&part;<sub>b</sub>&part;<sub>c</sub>h<sub>ad</sub> + &part;<sub>a</sub>&part;<sub>d</sub>h<sub>bc</sub>
&minus; &part;<sub>a</sub>&part;<sub>c</sub>h<sub>bd</sub> &minus; &part;<sub>b</sub>&part;<sub>d</sub>h<sub>ac</sub>)</code>
and indices are raised with <code>&eta;</code>
(<a href=\"modelica://GeneralRelativity.Perturbation.linearizedRiemann\">linearizedRiemann</a>,
<a href=\"modelica://GeneralRelativity.Perturbation.linearizedEinstein\">linearizedEinstein</a>).
In the Lorenz gauge the trace-reversed perturbation
(<a href=\"modelica://GeneralRelativity.Perturbation.traceReverse\">traceReverse</a>)
obeys <code>&#9633; h&#772;<sub>ab</sub> = -(16&pi;G/c<sup>4</sup>) T<sub>ab</sub></code>, which
<a href=\"modelica://GeneralRelativity.Numerical.LinearizedWave1D\">Numerical.LinearizedWave1D</a>
integrates.</p>
<h4>Newtonian limit</h4>
<p>For slow sources <code>h<sub>00</sub> = h<sub>ii</sub> = -2&Phi;/c&sup2;</code>
(<a href=\"modelica://GeneralRelativity.Perturbation.weakFieldMetric\">weakFieldMetric</a>)
and the <code>00</code> field equation is Poisson's equation
<code>&nabla;&sup2;&Phi; = 4&pi;G&rho;</code>. The tidal tensor
<code>R<sub>0i0j</sub> = &Phi;<sub>,ij</sub>/c&sup2;</code> comes from the Hessian of the potential
(<a href=\"modelica://GeneralRelativity.Perturbation.weakFieldCurvature\">weakFieldCurvature</a>).</p>
<h4>Why not difference the metric</h4>
<p>Earth's potential makes <code>h &asymp; 1.4&times;10<sup>-9</sup></code>. A double holds
sixteen digits, so second differences of <code>g</code> keep at most six or seven
digits of the curvature and first differences of the geodesic acceleration lose
three. The potentials in
<a href=\"modelica://GeneralRelativity.Perturbation.Potentials\">Potentials</a> are
analytic; their gradients and Hessians carry full precision, which is why the
Earth and cislunar examples integrate
<a href=\"modelica://GeneralRelativity.Perturbation.postNewtonianAcceleration\">postNewtonianAcceleration</a>
rather than the full geodesic equation with differenced Christoffel symbols.</p>
</html>"));
end Perturbation;
