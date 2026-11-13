within GeneralRelativity;
package NewtonCartan "Newtonian gravity as the geometry of an affine connection: clock form, degenerate spatial metric and the Newton-Cartan connection"
  annotation(Documentation(info="<html>
<p>Newton-Cartan theory writes Newtonian gravity in the language of general
relativity. Spacetime carries</p>
<ul>
<li>a <b>clock form</b> <code>&tau;<sub>a</sub> = &part;<sub>a</sub>(ct)</code>
(<a href=\"modelica://GeneralRelativity.NewtonCartan.clockForm\">clockForm</a>), the
absolute time;</li>
<li>a degenerate <b>spatial cometric</b> <code>h<sup>ab</sup></code> of rank 3 with
<code>h<sup>ab</sup>&tau;<sub>b</sub> = 0</code>
(<a href=\"modelica://GeneralRelativity.NewtonCartan.spatialCometric\">spatialCometric</a>), the
Euclidean geometry of each instant;</li>
<li>a torsion-free <b>affine connection</b> compatible with both,
<code>&nabla;&tau; = 0</code> and <code>&nabla;h = 0</code>. In an inertial or
rotating chart its only gravitational coefficient is
<code>&Gamma;<sup>i</sup><sub>00</sub> = &part;<sub>i</sub>&Phi;/c&sup2;</code>.</li>
</ul>
<p>Free fall is the geodesic equation of this connection, exactly as in general
relativity: with <code>&lambda; = ct</code>,
<code>d&sup2;x<sup>i</sup>/d&lambda;&sup2; = &minus;&Gamma;<sup>i</sup><sub>00</sub></code> is
<code>a = &minus;&nabla;&Phi;</code>. Gravity is curvature, not force: the Riemann tensor of
the connection is the tidal tensor, <code>R<sup>i</sup><sub>0j0</sub> = &part;<sub>i</sub>&part;<sub>j</sub>&Phi;/c&sup2;</code>,
and the field equation is <code>R<sub>ab</sub> = (4&pi;G/c&sup2;)&rho; &tau;<sub>a</sub>&tau;<sub>b</sub></code>
(<a href=\"modelica://GeneralRelativity.NewtonCartan.fieldEquationResidual\">fieldEquationResidual</a>),
which is Poisson's equation. Coriolis and centrifugal terms of a rotating
chart enter as <code>&Gamma;<sup>i</sup><sub>0j</sub></code> and leave the curvature zero.</p>
<h4>Relation to Schwarzschild</h4>
<p>In Schwarzschild coordinates <code>{ct, r, &theta;, &phi;}</code> the Levi-Civita
connection has <code>&Gamma;<sup>r</sup><sub>00</sub> = r<sub>s</sub>(1 &minus; r<sub>s</sub>/r)/(2r&sup2;)</code>,
while <a href=\"modelica://GeneralRelativity.NewtonCartan.pointMassSpherical\">pointMassSpherical</a>
has <code>r<sub>s</sub>/(2r&sup2;) = GM/(c&sup2;r&sup2;)</code> and keeps only the flat spatial
coefficients. The two agree to first order in <code>r<sub>s</sub>/r</code>
(the <code>c &rarr; &infin;</code> limit), so integrating both geodesic equations from
the same initial data isolates the relativistic corrections, as
<a href=\"modelica://GeneralRelativity.Examples.SchwarzschildNewtonCartanOrbit\">Examples.SchwarzschildNewtonCartanOrbit</a>
does. Geodesics of any connection are integrated by
<a href=\"modelica://GeneralRelativity.Geodesics.AffineGeodesic\">Geodesics.AffineGeodesic</a>;
curvature and torsion come from
<a href=\"modelica://GeneralRelativity.Connections.curvatureAt\">Connections.curvatureAt</a>.</p>
<p>Reference: E. Cartan (1923); K. Friedrichs (1927); Misner, Thorne and Wheeler,
<i>Gravitation</i>, chapter 12.</p>
</html>"));
end NewtonCartan;
