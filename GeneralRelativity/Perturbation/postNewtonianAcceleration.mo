within GeneralRelativity.Perturbation;
function postNewtonianAcceleration "First post-Newtonian correction to the acceleration of a test body about a point mass (harmonic coordinates, beta = gamma = 1)"
  input Real x[3](each unit = "m") "Position relative to the mass";
  input Real v[3](each unit = "m/s") "Coordinate velocity";
  input Real GM(unit = "m3/s2") "Gravitational parameter";
  input Real c(unit = "m/s") "Speed of light";
  output Real a[3](each unit = "m/s2") "Relativistic correction to add to the Newtonian acceleration";
protected
  Real r = sqrt(x * x);
algorithm
  a := GM / (c^2 * r^3) * ((4.0 * GM / r - v * v) * x + 4.0 * (x * v) * v);
  annotation(Documentation(info="<html>
<p>IERS Conventions (2010) equation 10.12 for a single central mass, the
Schwarzschild term of the Einstein-Infeld-Hoffmann equations. Adding it to
<code>-GM x/r&sup3;</code> reproduces perihelion advance of
<code>6&pi;GM/(c&sup2;a(1-e&sup2;))</code> per orbit. For Earth satellites it is of order
<code>10<sup>-9</sup></code> of the Newtonian acceleration, comparable to the
smallest geopotential terms used in precise orbit determination.</p>
</html>"));
end postNewtonianAcceleration;
