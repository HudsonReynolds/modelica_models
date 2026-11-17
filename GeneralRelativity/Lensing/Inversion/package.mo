within GeneralRelativity.Lensing;
package Inversion "Lens inversion: estimate the mass that distorts a lensed image, simulate the image it would make and refine the lens until the simulation matches"
  annotation(Documentation(info="<html>
<p>The forward problem (<a href=\"modelica://GeneralRelativity.Examples.GalaxyLensing\">Examples.GalaxyLensing</a>)
takes a lens mass and draws the lensed sky. This package runs it backwards. An image is a matrix of surface
brightness <code>image[i, j]</code> sampled at angles <code>(u[i], v[j])</code>, as in
<a href=\"modelica://GeneralRelativity.Visualization.writeFieldSlice\">writeFieldSlice</a>; the lens galaxy's own light
is assumed subtracted, as in real lens modelling. The steps:</p>
<ol>
<li><b>Input.</b> A matrix from a file (<code>Modelica.Utilities.Streams.readRealMatrix</code>; <code>tools/lens_image_to_mat.py</code>
converts a PNG or JPEG) or a synthetic image from <a href=\"modelica://GeneralRelativity.Lensing.Inversion.syntheticImage\">syntheticImage</a>.</li>
<li><b>Mass from the radius.</b> <a href=\"modelica://GeneralRelativity.Lensing.Inversion.einsteinRadiusFromImage\">einsteinRadiusFromImage</a>
fits a circle through the arcs; its radius is the Einstein radius and
<a href=\"modelica://GeneralRelativity.Lensing.massInEinsteinRadius\">massInEinsteinRadius</a> turns it into the projected mass inside
the ring, which does not depend on the mass profile.</li>
<li><b>Simulation.</b> <a href=\"modelica://GeneralRelativity.Lensing.Inversion.renderImage\">renderImage</a> draws the image that
a lens makes of S&eacute;rsic sources through the thin-lens equation (a point mass or an isothermal ellipse with external shear,
the deflections of <a href=\"modelica://GeneralRelativity.Lensing\">Lensing</a>) with sub-pixel sampling and a Gaussian
point-spread function. The simulation is the lightweight encoding of the picture: 7 lens numbers and 7 per source.</li>
<li><b>Comparison.</b> <a href=\"modelica://GeneralRelativity.Lensing.Inversion.imageResidual\">imageResidual</a>: sum of squares,
&chi;&sup2; for the pixel noise, the fraction of the image explained and the correlation.</li>
<li><b>Iteration.</b> <a href=\"modelica://GeneralRelativity.Lensing.Inversion.fitLens\">fitLens</a> refines lens and source by
Levenberg-Marquardt least squares on the pixels until the simulated image matches. Gaussian priors on any parameter turn the
step into the iterated extended Kalman measurement update (a maximum a-posteriori estimate): this is where outside knowledge enters,
for example the position or mass of a catalogued galaxy. The fit returns the posterior covariance, and
<a href=\"modelica://GeneralRelativity.Lensing.Inversion.einsteinMass\">einsteinMass</a> turns the fitted lens into the mass inside
its critical curve with a standard error.</li>
</ol>
<p><a href=\"modelica://GeneralRelativity.Lensing.Inversion.invert\">invert</a> runs steps 2 to 5 with several starting points
(pixel fits of lenses have local minima: an off-axis source, an unknown orientation, a compact source), and
<a href=\"modelica://GeneralRelativity.Visualization.writeLensInversion\">Visualization.writeLensInversion</a> runs it and writes the
result for the viewer.</p>
<h4>Parameter vectors</h4>
<table>
<tr><th>Lens <code>lens[7]</code></th><th>Source row <code>sources[k, 7]</code></th></tr>
<tr><td>1 deflection scale b (rad): the semi-major axis of the critical curve</td><td>1 x (rad)</td></tr>
<tr><td>2 x of the lens centre (rad)</td><td>2 y (rad)</td></tr>
<tr><td>3 y of the lens centre (rad)</td><td>3 brightness at the half-light radius I<sub>e</sub></td></tr>
<tr><td>4 axis ratio q (ignored by a point mass)</td><td>4 half-light radius r<sub>e</sub> (rad)</td></tr>
<tr><td>5 major-axis position angle (rad)</td><td>5 axis ratio q</td></tr>
<tr><td>6 external shear &gamma;</td><td>6 position angle (rad)</td></tr>
<tr><td>7 shear position angle (rad)</td><td>7 S&eacute;rsic index n (1 exponential disk, 4 de Vaucouleurs)</td></tr>
</table>
<h4>Which radius, which mass</h4>
<p>For the isothermal ellipse the critical curve has semi-axes <code>b</code> and <code>qb</code>. Its equal-area radius
<code>&theta;<sub>E</sub> = b&radic;q</code> is what a circle fitted to the arcs measures, and the projected mass inside it is
<code>&pi;qb&sup2;&Sigma;<sub>cr</sub>D<sub>l</sub>&sup2;</code> = <code>massInEinsteinRadius(&theta;<sub>E</sub>)</code>
(<a href=\"modelica://GeneralRelativity.Lensing.Inversion.einsteinMass\">einsteinMass</a>). For a point mass both are just
<code>b</code> and the total mass. The radius-based estimate of step 2 and the fit therefore report the same quantity.</p>
<h4>Examples</h4>
<ul>
<li><a href=\"modelica://GeneralRelativity.Examples.LensInversion\">Examples.LensInversion</a>: one synthetic image whose true mass
is known; the recovery is asserted.</li>
<li><a href=\"modelica://GeneralRelativity.Examples.LensInversionLibrary\">Examples.LensInversionLibrary</a>: eight synthetic images
that look different on the sky (ring, double, giant arc, Einstein cross, point mass, noisy, group-scale, low-mass), one scene with a
gallery to choose from.</li>
<li><a href=\"modelica://GeneralRelativity.Examples.LensInversionFromFile\">Examples.LensInversionFromFile</a>: an image read from a
file, true lens unknown.</li>
</ul>
<p>Each writes a scene for the viewer's <b>lensfit</b> style (observed, first and fitted simulations, residual, &chi;&sup2; and
mass per iteration, parameter table) and a <code>lens</code> record so the <b>lensing</b> style can show a sky through the recovered
lens.</p>
<p>Limitations: one lens plane, parametric sources (S&eacute;rsic profiles), a Gaussian point-spread function, no pixelated source
reconstruction, no lens-light model. Angles are radians; multiply by <code>Constants.arcsec</code> to pass arcseconds.</p>
</html>"));
end Inversion;
