within GeneralRelativity.Numerical;
model LinearizedWave1D "Method-of-lines wave equation for a trace-reversed metric perturbation component with outflow boundaries"
  parameter Integer n = 201 "Number of nodes";
  parameter Real length = 100.0 "Domain length, x in [0, length]";
  parameter Real c = 1.0 "Wave speed";
  parameter Real amplitude = 1.0e-3 "Initial pulse amplitude";
  parameter Real centre = 50.0 "Initial pulse centre";
  parameter Real width = 4.0 "Initial pulse standard deviation";
  parameter Real sourceAmplitude = 0.0 "Amplitude of a stationary Gaussian source -16 pi G/c^4 T";
  parameter Real sourceCentre = 50.0 "Source centre";
  parameter Real sourceWidth = 4.0 "Source standard deviation";
  final parameter Real dx = length / (n - 1) "Node spacing";
  final parameter Real xNode[n] = {(i - 1) * dx for i in 1:n} "Node positions";
  final parameter Real hStart[n] = {amplitude * exp(-0.5 * ((xNode[i] - centre) / width)^2) for i in 1:n} "Initial pulse";
  final parameter Real initialEnergy = 0.5 * dx * c^2 * sum(((hStart[i + 1] - hStart[i]) / dx)^2 for i in 1:n - 1)
    "Discrete energy of the initial pulse";
  Real h[n](start = hStart, each fixed = true) "Perturbation hbar at the nodes";
  Real v[n](each start = 0.0) "d h / d t";
  Real energy "Discrete wave energy sum(v^2 + c^2 (dh/dx)^2) dx / 2, non-increasing with outflow boundaries";
  Real hMax "Largest |h| on the grid";
initial equation
  for i in 2:n - 1 loop
    v[i] = 0.0;
  end for;
equation
  for i in 2:n - 1 loop
    der(h[i]) = v[i];
    der(v[i]) = c^2 * (h[i + 1] - 2.0 * h[i] + h[i - 1]) / dx^2
      + sourceAmplitude * exp(-0.5 * ((xNode[i] - sourceCentre) / sourceWidth)^2);
  end for;
  // Sommerfeld radiation conditions: left-going at x = 0, right-going at x = length.
  der(h[1]) = c * (h[2] - h[1]) / dx;
  der(h[n]) = -c * (h[n] - h[n - 1]) / dx;
  v[1] = der(h[1]);
  v[n] = der(h[n]);
  energy = 0.5 * dx * (sum(v[i]^2 for i in 1:n) + c^2 * sum(((h[i + 1] - h[i]) / dx)^2 for i in 1:n - 1));
  hMax = max(abs(h));
  annotation(Documentation(info="<html>
<p>In the Lorenz gauge each component of the trace-reversed perturbation obeys
<code>&#9633; h&#772; = -(16&pi;G/c<sup>4</sup>) T</code>. This model keeps one
component on a line with second-order central differences in space and leaves
time to the Modelica integrator. A Gaussian pulse splits into two halves that
travel at <code>c</code> and leave through the boundaries; a stationary source
builds up the static field <code>h&#772; &prop; -&Phi;</code>. Use an interval of
about <code>dx/c</code> for output and a solver tolerance of <code>1e-8</code>.</p>
</html>"));
end LinearizedWave1D;
