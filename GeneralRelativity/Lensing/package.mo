within GeneralRelativity;
package Lensing "Gravitational lensing: cosmological distances, Einstein radii and thin-lens deflections of point masses and galaxies"
  annotation(Documentation(info="<html>
<p>Thin-lens gravitational lensing. Light from a source at angular position <b>&beta;</b> reaches the
observer from the image positions <b>&theta;</b> that solve the lens equation</p>
<p><code>&beta; = &theta; &minus; &alpha;(&theta;)</code>,</p>
<p>where the scaled deflection <b>&alpha;</b> is the true bending angle times <code>D<sub>ls</sub>/D<sub>s</sub></code>
(distances to the lens, to the source and between them, all angular-diameter distances).
For a point mass the bending angle is <code>4GM/(c&sup2;b)</code>, the weak-field limit of a
Schwarzschild null geodesic
(<a href=\"modelica://GeneralRelativity.Examples.LightDeflection\">Examples.LightDeflection</a> integrates one and
checks it). Angles are in radians unless a name says otherwise.</p>
<table>
<tr><th>Function</th><th>Gives</th></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Lensing.angularDiameterDistance\">angularDiameterDistance</a></td><td><code>D<sub>A</sub>(z<sub>1</sub>, z<sub>2</sub>)</code> in a Friedmann-Lema&icirc;tre universe with matter and &Lambda;</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Lensing.criticalDensity\">criticalDensity</a></td><td><code>&Sigma;<sub>cr</sub> = c&sup2;D<sub>s</sub>/(4&pi;G D<sub>l</sub>D<sub>ls</sub>)</code>; the convergence is <code>&kappa; = &Sigma;/&Sigma;<sub>cr</sub></code></td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Lensing.einsteinRadiusPointMass\">einsteinRadiusPointMass</a>, <a href=\"modelica://GeneralRelativity.Lensing.einsteinRadiusSIS\">einsteinRadiusSIS</a></td><td>Einstein radii</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Lensing.deflectionPointMass\">deflectionPointMass</a>, <a href=\"modelica://GeneralRelativity.Lensing.deflectionSIE\">deflectionSIE</a>, <a href=\"modelica://GeneralRelativity.Lensing.deflectionShear\">deflectionShear</a></td><td>Scaled deflections: point mass, singular isothermal ellipse (a galaxy), external shear</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Lensing.convergenceSIE\">convergenceSIE</a></td><td>&kappa; of the isothermal ellipse</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Lensing.pointMassImages\">pointMassImages</a></td><td>Both images of a point-mass lens and their magnifications</td></tr>
</table>
<p><a href=\"modelica://GeneralRelativity.Examples.GalaxyLensing\">Examples.GalaxyLensing</a> writes a scene for the
viewer's <b>lensing</b> style: a telescope-like image of lensed background galaxies whose lens mass,
shape and redshifts are sliders.</p>
</html>"));
end Lensing;
