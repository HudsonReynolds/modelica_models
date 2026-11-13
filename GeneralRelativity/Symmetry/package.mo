within GeneralRelativity;
package Symmetry "Lie derivative of the metric, Killing and conformal Killing tests, and standard symmetry generators"
  annotation(Documentation(info="<html>
<p>A vector field <code>&xi;</code> is a <b>Killing vector</b> when the metric is
invariant along its flow, <code>&Lscr;<sub>&xi;</sub>g = 0</code>. In coordinates</p>
<pre>(&Lscr;<sub>&xi;</sub>g)<sub>ab</sub> = &xi;<sup>c</sup>&part;<sub>c</sub>g<sub>ab</sub> + g<sub>cb</sub>&part;<sub>a</sub>&xi;<sup>c</sup> + g<sub>ac</sub>&part;<sub>b</sub>&xi;<sup>c</sup>
                = &nabla;<sub>a</sub>&xi;<sub>b</sub> + &nabla;<sub>b</sub>&xi;<sub>a</sub>.</pre>
<p>A <b>conformal Killing vector</b> satisfies
<code>&Lscr;<sub>&xi;</sub>g = (2/n)(&nabla;<sub>c</sub>&xi;<sup>c</sup>) g</code> with
<code>n = 4</code>. Every Killing vector gives a constant of geodesic motion,
<code>g(&xi;, u)</code>
(<a href=\"modelica://GeneralRelativity.Geodesics.killingConserved\">Geodesics.killingConserved</a>).</p>
<p>Vector fields have the signature of
<a href=\"modelica://GeneralRelativity.Interfaces.VectorFieldFunction\">Interfaces.VectorFieldFunction</a>,
<code>(x[4], p[:]) &rarr; &xi;[4]</code>, and are passed by name like metrics. The operators take the metric parameters as <code>p</code> and the field parameters as <code>q</code>; the table lists the latter. Their
derivatives and the metric derivatives are central finite differences with a
common step per coordinate.
<a href=\"modelica://GeneralRelativity.Symmetry.killingResidual\">killingResidual</a>
divides <code>max|&Lscr;<sub>&xi;</sub>g|</code> by the size of the individual
terms, so a Killing vector gives a value of the order of the finite-difference
error (about 1e-7 for step 1e-3) and a generic field a value of order one.</p>
<h4>Generators</h4>
<table>
<tr><th>Field</th><th>Chart</th><th>Components</th></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Symmetry.coordinateVector\">coordinateVector</a></td><td>any</td><td><code>&part;<sub>k</sub></code>, <code>q = {k}</code></td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Symmetry.timeTranslation\">timeTranslation</a></td><td>any</td><td><code>&part;<sub>0</sub></code></td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Symmetry.axialRotation\">axialRotation</a></td><td>spherical {ct, r, &theta;, &phi;}</td><td><code>&part;<sub>&phi;</sub></code></td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Symmetry.rotationX\">rotationX</a>, <a href=\"modelica://GeneralRelativity.Symmetry.rotationY\">rotationY</a></td><td>spherical</td><td>rotations about the x and y axes</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Symmetry.cartesianTranslation\">cartesianTranslation</a></td><td>Cartesian {ct, x, y, z}</td><td>constant <code>q</code></td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Symmetry.cartesianRotation\">cartesianRotation</a></td><td>Cartesian</td><td><code>q &times; r</code></td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Symmetry.cartesianBoost\">cartesianBoost</a></td><td>Cartesian</td><td><code>(q&middot;r) &part;<sub>0</sub> + x<sup>0</sup> q</code></td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Symmetry.dilation\">dilation</a></td><td>Cartesian</td><td><code>x<sup>a</sup>&part;<sub>a</sub></code>, conformal Killing in Minkowski</td></tr>
</table>
<p>Schwarzschild has <code>&part;<sub>t</sub></code> and the three rotations; Kerr
keeps only <code>&part;<sub>t</sub></code> and <code>&part;<sub>&phi;</sub></code>;
Minkowski has the ten Poincar&eacute; generators.</p>
</html>"));
end Symmetry;
