within GeneralRelativity.Lensing;
package Inversion "Lens inversion: estimate the mass that distorts a lensed image, simulate the image it would make and refine the lens until the simulation matches"
  annotation(Documentation(info="<html>
<p>The forward problem (<a href=\"modelica://GeneralRelativity.Examples.GalaxyLensing\">Examples.GalaxyLensing</a>)
takes a lens mass and draws the lensed sky. This package runs it backwards. An image is a matrix of surface
brightness <code>image[i, j]</code> sampled at angles <code>(u[i], v[j])</code>, as in
<a href=\"modelica://GeneralRelativity.Visualization.writeFieldSlice\">writeFieldSlice</a>; the lens galaxy's own light
is assumed subtracted, as in real lens modelling. The steps, each a function here:</p>
<ol>
<li><b>Input.</b> A matrix from a file (<code>Modelica.Utilities.Streams.readRealMatrix</code>; <code>tools/lens_image_to_mat.py</code>
converts a PNG or JPEG) or a synthetic image from <a href=\"modelica://GeneralRelativity.Lensing.Inversion.syntheticImage\">syntheticImage</a>.</li>
<li><b>Mass from the radius.</b> <a href=\"modelica://GeneralRelativity.Lensing.Inversion.brightnessCentroid\">brightnessCentroid</a> locates
the lens and <a href=\"modelica://GeneralRelativity.Lensing.Inversion.einsteinRadiusFromImage\">einsteinRadiusFromImage</a> reads the ring
or arc radius; <a href=\"modelica://GeneralRelativity.Lensing.massInEinsteinRadius\">massInEinsteinRadius</a> turns it into the projected
mass inside the ring, which is model-independent, and
<a href=\"modelica://GeneralRelativity.Lensing.velocityDispersionSIS\">velocityDispersionSIS</a> into a galaxy's velocity dispersion.</li>
<li><b>Simulation.</b> <a href=\"modelica://GeneralRelativity.Lensing.Inversion.renderImage\">renderImage</a> draws the image that
lens makes of a S&eacute;rsic source through the thin-lens equation (a point mass or an isothermal ellipse with external shear,
the deflections of <a href=\"modelica://GeneralRelativity.Lensing\">Lensing</a>) with sub-pixel sampling and an optional Gaussian
point-spread function. The simulation is the lightweight encoding of the picture: 7 lens numbers and 7 per source.</li>
<li><b>Comparison.</b> <a href=\"modelica://GeneralRelativity.Lensing.Inversion.imageResidual\">imageResidual</a>: sum of squares,
&chi;&sup2; for the pixel noise, the fraction of the image explained and the correlation.</li>
<li><b>Iteration.</b> <a href=\"modelica://GeneralRelativity.Lensing.Inversion.fitLens\">fitLens</a> refines lens and source by
Levenberg-Marquardt least squares on the pixels until the simulated image matches. Gaussian priors on any parameter turn the
step into the iterated extended Kalman measurement update (a maximum a-posteriori estimate): this is where outside knowledge enters,
for example the position or mass of a catalogued galaxy. The fit also returns the standard error of every parameter, hence of the mass.</li>
</ol>
<h4>Parameter vectors</h4>
<table>
<tr><th>Lens <code>lens[7]</code></th><th>Source row <code>sources[k, 7]</code></th></tr>
<tr><td>1 Einstein radius &theta;<sub>E</sub> (rad)</td><td>1 x (rad)</td></tr>
<tr><td>2 x of the lens centre (rad)</td><td>2 y (rad)</td></tr>
<tr><td>3 y of the lens centre (rad)</td><td>3 brightness at the half-light radius I<sub>e</sub></td></tr>
<tr><td>4 axis ratio q (ignored by a point mass)</td><td>4 half-light radius r<sub>e</sub> (rad)</td></tr>
<tr><td>5 major-axis position angle (rad)</td><td>5 axis ratio q</td></tr>
<tr><td>6 external shear &gamma;</td><td>6 position angle (rad)</td></tr>
<tr><td>7 shear position angle (rad)</td><td>7 S&eacute;rsic index n (1 exponential disk, 4 de Vaucouleurs)</td></tr>
</table>
<p><a href=\"modelica://GeneralRelativity.Examples.LensInversion\">Examples.LensInversion</a> runs all five steps on a synthetic
image whose true mass it knows and asserts the recovery;
<a href=\"modelica://GeneralRelativity.Examples.LensInversionFromFile\">Examples.LensInversionFromFile</a> does the same for a
matrix read from a file. Both write a scene for the viewer's <b>lensfit</b> style (observed, simulated and residual images, the
&chi;&sup2; and mass history) and a <code>lens</code> record so the <b>lensing</b> style can show a sky through the recovered lens.</p>
<p>Limitations: one lens plane, parametric sources (several S&eacute;rsic profiles), a Gaussian point-spread function, no
pixelated source reconstruction, no lens-light model. Angles are radians; multiply by
<code>Constants.arcsec</code> to pass arcseconds.</p>
</html>"));
end Inversion;
