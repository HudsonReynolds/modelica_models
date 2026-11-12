within GeneralRelativity.Geodesics;
model Geodesic "Geodesic of a replaceable metric, integrated in the affine parameter (Modelica time)"
  replaceable function metric = GeneralRelativity.Metrics.minkowski
    constrainedby GeneralRelativity.Interfaces.MetricFunction "Metric g_ij(x, p)";
  parameter Real p[:] = {0.0} "Metric parameter vector in the layout the metric documents";
  parameter Real c = 1.0 "Speed of light in the chart's units; 1 for geometric units";
  parameter Real x_start[4] = {0.0, 0.0, 0.0, 0.0} "Initial coordinates; x[1] = c*t";
  parameter Real u_start[4] = {1.0, 0.0, 0.0, 0.0} "Initial four-velocity dx/d lambda";
  parameter Real step[4] = {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3} "Finite-difference step per coordinate for metric gradients";
  parameter Boolean massive = true "true: lambda is proper time and g(u,u) = -c^2; false: null geodesic";
  Real x[4](start = x_start, each fixed = true) "Coordinates";
  Real u[4](start = u_start, each fixed = true) "Four-velocity dx/d lambda";
  Real a[4] "Geodesic acceleration";
  Real g[4, 4] "Metric at the particle";
  Real Gamma[4, 4, 4] "Christoffel symbols at the particle";
  Real normResidual "g(u,u) + c^2 (massive) or g(u,u) (null); drift measures integration error";
  Real t "Coordinate time x[1]/c";
equation
  (a, g, Gamma) = GeneralRelativity.Geodesics.geodesicAcceleration(metric, x, p, step, u);
  der(x) = u;
  der(u) = a;
  normResidual = GeneralRelativity.Geodesics.normResidual(g, u, if massive then c else 0.0);
  t = x[1] / c;
  annotation(Documentation(info="<html>
<p>Redeclare the metric and bind its parameters:</p>
<pre>model Orbit
  extends GeneralRelativity.Geodesics.Geodesic(
    redeclare function metric = GeneralRelativity.Metrics.schwarzschild,
    p = {2.0}, x_start = {0, 10, Modelica.Constants.pi/2, 0},
    u_start = {E/(1 - 2/10), 0, 0, L/100});
end Orbit;</pre>
<p>Modelica <code>time</code> is the affine parameter: proper time for a massive
particle. Use a tight tolerance (<code>1e-10</code>) for long integrations and
watch <code>normResidual</code>, which an exact integrator keeps at zero.
The Christoffel symbols come from central differences of the metric with
<code>step</code>; choose it about <code>10<sup>-4</sup></code> of each coordinate's
scale and smaller than the distance to any chart singularity.</p>
</html>"));
end Geodesic;
