within GeneralRelativity.Geodesics;
model AffineGeodesic "Autoparallel of a replaceable affine connection, such as the Newton-Cartan connection, with a selectable clock rate"
  replaceable function connection = GeneralRelativity.NewtonCartan.pointMassCartesian
    constrainedby GeneralRelativity.Interfaces.ConnectionFunction "Connection Gamma^i_jk(x, p)";
  parameter Real p[:] = {0.0} "Connection parameter vector in the layout the connection documents";
  parameter Real c = 1.0 "Speed of light in the chart's units, for t = x[1]/c";
  parameter Real x_start[4] = {0.0, 0.0, 0.0, 0.0} "Initial coordinates; x[1] = c*t";
  parameter Real u_start[4] = {1.0, 0.0, 0.0, 0.0} "Initial tangent dx/d lambda";
  Real lambdaRate = 1.0 "d lambda / d time: 1 integrates in the affine parameter itself; bind it to synchronise with another clock";
  Real x[4](start = x_start, each fixed = true) "Coordinates";
  Real u[4](start = u_start, each fixed = true) "Tangent dx/d lambda";
  Real a[4] "d^2 x / d lambda^2 = -Gamma^i_jk u^j u^k";
  Real Gamma[4, 4, 4] "Connection coefficients at the particle";
  Real t "Coordinate time x[1]/c";
equation
  Gamma = connection(x, p);
  a = GeneralRelativity.Geodesics.acceleration(Gamma, u);
  der(x) = lambdaRate * u;
  der(u) = lambdaRate * a;
  t = x[1] / c;
  annotation(Documentation(info="<html>
<p>Integrates <code>d&sup2;x<sup>i</sup>/d&lambda;&sup2; = &minus;&Gamma;<sup>i</sup><sub>jk</sub>u<sup>j</sup>u<sup>k</sup></code>
for any <a href=\"modelica://GeneralRelativity.Interfaces.ConnectionFunction\">ConnectionFunction</a>;
only the symmetric part of <code>&Gamma;</code> enters, so with torsion this is the
autoparallel. Unlike <a href=\"modelica://GeneralRelativity.Geodesics.Geodesic\">Geodesic</a>
there is no metric, hence no normalization of <code>u</code>.</p>
<p>Modelica <code>time</code> need not be <code>&lambda;</code>: the model solves
<code>dx/d(time) = lambdaRate&middot;u</code>. For a Newton-Cartan body take
<code>&lambda; = ct</code> (so <code>u<sup>0</sup> = 1</code>, and <code>u<sup>0</sup></code> stays 1
because <code>&Gamma;<sup>0</sup><sub>ab</sub> = 0</code>) and bind
<code>lambdaRate</code> to <code>c&nbsp;dt/d(time)</code> of another body to advance both
in the same coordinate time; see
<a href=\"modelica://GeneralRelativity.Examples.SchwarzschildNewtonCartanOrbit\">Examples.SchwarzschildNewtonCartanOrbit</a>.</p>
</html>"));
end AffineGeodesic;
