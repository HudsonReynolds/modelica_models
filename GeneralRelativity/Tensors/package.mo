within GeneralRelativity;
package Tensors "Index algebra: inverse metric, Christoffel symbols, Riemann, Ricci, Einstein and stress-energy tensors"
  annotation(Documentation(info="<html>
<p>These functions are pure algebra on a metric and its coordinate derivatives,
so they accept analytic derivatives as well as the finite differences produced
by <a href=\"modelica://GeneralRelativity.Metrics.derivative\">Metrics.derivative</a>
and <a href=\"modelica://GeneralRelativity.Metrics.secondDerivative\">Metrics.secondDerivative</a>.
Call them in order:</p>
<ol>
<li><a href=\"modelica://GeneralRelativity.Tensors.inverseMetric\">inverseMetric</a>: <code>g<sup>ij</sup></code>;</li>
<li><a href=\"modelica://GeneralRelativity.Tensors.christoffel\">christoffel</a>: <code>&Gamma;<sup>i</sup><sub>jk</sub></code> from <code>g</code>, <code>dg</code>;</li>
<li><a href=\"modelica://GeneralRelativity.Tensors.riemannLower\">riemannLower</a>: <code>R<sub>abcd</sub></code> from <code>g</code>, <code>ddg</code>, <code>&Gamma;</code>;</li>
<li><a href=\"modelica://GeneralRelativity.Tensors.ricci\">ricci</a>, <a href=\"modelica://GeneralRelativity.Tensors.ricciScalar\">ricciScalar</a>,
<a href=\"modelica://GeneralRelativity.Tensors.einstein\">einstein</a>,
<a href=\"modelica://GeneralRelativity.Tensors.stressEnergyFromEinstein\">stressEnergyFromEinstein</a>;</li>
<li>invariants: <a href=\"modelica://GeneralRelativity.Tensors.kretschmann\">kretschmann</a>.</li>
</ol>
<p><a href=\"modelica://GeneralRelativity.Tensors.Curvature\">Curvature</a> bundles
all of them for one point;
<a href=\"modelica://GeneralRelativity.Metrics.curvatureAt\">Metrics.curvatureAt</a>
fills it from a metric function.</p>
<p>The Riemann tensor is built in its fully covariant form
<code>R<sub>abcd</sub> = &frac12;(&part;<sub>b</sub>&part;<sub>c</sub>g<sub>ad</sub> + &part;<sub>a</sub>&part;<sub>d</sub>g<sub>bc</sub>
&minus; &part;<sub>a</sub>&part;<sub>c</sub>g<sub>bd</sub> &minus; &part;<sub>b</sub>&part;<sub>d</sub>g<sub>ac</sub>)
+ g<sub>ef</sub>(&Gamma;<sup>e</sup><sub>bc</sub>&Gamma;<sup>f</sup><sub>ad</sub> &minus; &Gamma;<sup>e</sup><sub>bd</sub>&Gamma;<sup>f</sup><sub>ac</sub>)</code>,
which needs only one layer of differentiation of the metric and has the pair
symmetries exactly. Its sign agrees with the mixed convention in the package
documentation: on the unit sphere <code>R<sub>&theta;&phi;&theta;&phi;</sub> = sin&sup2;&theta;</code>.</p>
</html>"));
end Tensors;
