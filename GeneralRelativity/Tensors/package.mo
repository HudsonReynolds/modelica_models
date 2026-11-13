within GeneralRelativity;
package Tensors "Index algebra: inverse metric, Christoffel symbols, Riemann, Ricci, Einstein and stress-energy tensors"
  annotation(Documentation(info="<html>
<p>These functions are pure algebra on a metric and its coordinate derivatives,
so they accept analytic derivatives as well as the finite differences produced
by <a href=\"modelica://GeneralRelativity.Metrics.derivative\">Metrics.derivative</a>
and <a href=\"modelica://GeneralRelativity.Metrics.secondDerivative\">Metrics.secondDerivative</a>.
The core chain below is dimension-generic: array sizes follow the metric, so
the same functions give the curvature of a 2-sphere, of the spatial 3-metric
of a 3+1 split, or of a spacetime. Functions built on the 4-dimensional
Levi-Civita symbol or Weyl decomposition say so. Call them in order:</p>
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
<p>Curvature decomposition and invariants:
<a href=\"modelica://GeneralRelativity.Tensors.schouten\">schouten</a> (<code>P<sub>ab</sub></code>),
<a href=\"modelica://GeneralRelativity.Tensors.weyl\">weyl</a> (<code>C<sub>abcd</sub> = R<sub>abcd</sub> &minus; g &#8743; P</code>) and
<a href=\"modelica://GeneralRelativity.Tensors.traceFreeRicci\">traceFreeRicci</a> split the Riemann tensor into
its irreducible parts in any dimension;
<a href=\"modelica://GeneralRelativity.Tensors.ricciSquared\">ricciSquared</a>,
<a href=\"modelica://GeneralRelativity.Tensors.fullContraction\">fullContraction</a> (for <code>C<sub>abcd</sub>C<sup>abcd</sup></code>),
<a href=\"modelica://GeneralRelativity.Tensors.gaussBonnet\">gaussBonnet</a>,
<a href=\"modelica://GeneralRelativity.Tensors.cubicInvariants\">cubicInvariants</a> and
<a href=\"modelica://GeneralRelativity.Tensors.weylCubed\">weylCubed</a> (built on
<a href=\"modelica://GeneralRelativity.Tensors.bivectorMatrix\">bivectorMatrix</a>) give the scalar invariants,
collected in an <a href=\"modelica://GeneralRelativity.Tensors.Invariants\">Invariants</a> record by
<a href=\"modelica://GeneralRelativity.Tensors.curvatureInvariants\">curvatureInvariants</a>.
The Cotton tensor needs derivatives of the Schouten tensor and lives in
<a href=\"modelica://GeneralRelativity.Metrics.cottonAt\">Metrics.cottonAt</a>; tetrads, Weyl scalars,
Petrov types and energy conditions are in <a href=\"modelica://GeneralRelativity.Frames\">Frames</a>.</p>
<p>The Riemann tensor is built in its fully covariant form
<code>R<sub>abcd</sub> = &frac12;(&part;<sub>b</sub>&part;<sub>c</sub>g<sub>ad</sub> + &part;<sub>a</sub>&part;<sub>d</sub>g<sub>bc</sub>
&minus; &part;<sub>a</sub>&part;<sub>c</sub>g<sub>bd</sub> &minus; &part;<sub>b</sub>&part;<sub>d</sub>g<sub>ac</sub>)
+ g<sub>ef</sub>(&Gamma;<sup>e</sup><sub>bc</sub>&Gamma;<sup>f</sup><sub>ad</sub> &minus; &Gamma;<sup>e</sup><sub>bd</sub>&Gamma;<sup>f</sup><sub>ac</sub>)</code>,
which needs only one layer of differentiation of the metric and has the pair
symmetries exactly. Its sign agrees with the mixed convention in the package
documentation: on the unit sphere <code>R<sub>&theta;&phi;&theta;&phi;</sub> = sin&sup2;&theta;</code>.</p>
</html>"));
end Tensors;
