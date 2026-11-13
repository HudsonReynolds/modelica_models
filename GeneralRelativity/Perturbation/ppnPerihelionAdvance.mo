within GeneralRelativity.Perturbation;
function ppnPerihelionAdvance "Relativistic perihelion advance per orbit, (2 + 2 gamma - beta)/3 * 6 pi GM / (c^2 a (1 - e^2))"
  input Real GM(unit = "m3/s2") "Gravitational parameter of the central body (or of the pair)";
  input Real a(unit = "m") "Semi-major axis";
  input Real e "Eccentricity";
  input Real gamma = 1.0 "PPN gamma";
  input Real beta = 1.0 "PPN beta";
  input Real c(unit = "m/s") = GeneralRelativity.Constants.c "Speed of light";
  output Real deltaOmega(unit = "rad") "Advance of the argument of perihelion per orbit";
algorithm
  deltaOmega := (2.0 + 2.0 * gamma - beta) / 3.0 * 6.0 * GeneralRelativity.Constants.pi * GM / (c^2 * a * (1.0 - e^2));
  annotation(Documentation(info="<html>
<p>The point-mass part of the PPN perihelion shift (Will, Living Rev. Relativ. 17:4, 2014) for a test
body; preferred-frame terms and the quadrupole moment of the central body are not
included. Multiply by the number of orbits per century and by
<code>180&times;3600/&pi;</code> for arcseconds per century: Mercury gives 42.98&Prime;/cy in
general relativity. See
<a href=\"modelica://GeneralRelativity.Examples.SchwarzschildOrbit\">Examples.SchwarzschildOrbit</a>
for the same effect from the full geodesic equation.</p>
</html>"));
end ppnPerihelionAdvance;
