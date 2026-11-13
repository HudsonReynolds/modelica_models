within GeneralRelativity;
package Connections "General affine connections: curvature and torsion of any connection, non-metricity, contortion and the Levi-Civita connection of a metric"
  annotation(Documentation(info="<html>
<p>A connection on the tangent bundle is given by its coefficients in a chart,
<code>&nabla;<sub>j</sub>v<sup>i</sup> = &part;<sub>j</sub>v<sup>i</sup> + &Gamma;<sup>i</sup><sub>jk</sub>v<sup>k</sup></code>,
stored as <code>Gamma[i, j, k]</code> with the <b>derivative index second</b>, as for the
Christoffel symbols of <a href=\"modelica://GeneralRelativity.Tensors.christoffel\">Tensors.christoffel</a>.
Nothing here assumes the connection comes from a metric, so the same functions
handle the Levi-Civita connection of general relativity, the Newton-Cartan
connection of Newtonian gravity
(<a href=\"modelica://GeneralRelativity.NewtonCartan\">NewtonCartan</a>) and
metric-affine connections with torsion and non-metricity.</p>
<p>A connection is passed as a function with the signature of
<a href=\"modelica://GeneralRelativity.Interfaces.ConnectionFunction\">Interfaces.ConnectionFunction</a>;
<a href=\"modelica://GeneralRelativity.Connections.derivative\">derivative</a> differences it, and
<a href=\"modelica://GeneralRelativity.Connections.curvatureAt\">curvatureAt</a> returns its
Riemann tensor, Ricci tensor and torsion at a point.</p>
<h4>Conventions</h4>
<ul>
<li>torsion <code>T<sup>i</sup><sub>jk</sub> = &Gamma;<sup>i</sup><sub>jk</sub> &minus; &Gamma;<sup>i</sup><sub>kj</sub></code>;</li>
<li>curvature <code>R<sup>&rho;</sup><sub>&sigma;&mu;&nu;</sub> = &part;<sub>&mu;</sub>&Gamma;<sup>&rho;</sup><sub>&nu;&sigma;</sub> &minus; &part;<sub>&nu;</sub>&Gamma;<sup>&rho;</sup><sub>&mu;&sigma;</sub> + &Gamma;<sup>&rho;</sup><sub>&mu;&lambda;</sub>&Gamma;<sup>&lambda;</sup><sub>&nu;&sigma;</sub> &minus; &Gamma;<sup>&rho;</sup><sub>&nu;&lambda;</sub>&Gamma;<sup>&lambda;</sup><sub>&mu;&sigma;</sub></code>,
so <code>[&nabla;<sub>&mu;</sub>, &nabla;<sub>&nu;</sub>]v<sup>&rho;</sup> = R<sup>&rho;</sup><sub>&sigma;&mu;&nu;</sub>v<sup>&sigma;</sup> &minus; T<sup>&lambda;</sup><sub>&mu;&nu;</sub>&nabla;<sub>&lambda;</sub>v<sup>&rho;</sup></code>;
for a Levi-Civita connection this is the package-wide Riemann convention;</li>
<li>Ricci <code>R<sub>&sigma;&nu;</sub> = R<sup>&rho;</sup><sub>&sigma;&rho;&nu;</sub></code>, not symmetric in general;</li>
<li>non-metricity <code>Q<sub>&lambda;&mu;&nu;</sub> = &nabla;<sub>&lambda;</sub>g<sub>&mu;&nu;</sub></code>;</li>
<li>contortion <code>K<sup>i</sup><sub>jk</sub> = &Gamma;<sup>i</sup><sub>jk</sub> &minus; {<sup>i</sup><sub>jk</sub>}</code>
for a metric-compatible connection; with indices lowered on the first slot
<code>K<sub>abc</sub> = &frac12;(T<sub>abc</sub> &minus; T<sub>bca</sub> + T<sub>cab</sub>)</code>.</li>
</ul>
</html>"));
end Connections;
