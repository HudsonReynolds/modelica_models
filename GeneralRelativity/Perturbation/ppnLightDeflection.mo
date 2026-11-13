within GeneralRelativity.Perturbation;
function ppnLightDeflection "Deflection angle of light passing a mass, (1 + gamma)/2 * 4 GM / (c^2 b)"
  input Real GM(unit = "m3/s2") "Gravitational parameter of the deflector";
  input Real b(unit = "m") "Impact parameter (closest approach)";
  input Real gamma = 1.0 "PPN gamma";
  input Real c(unit = "m/s") = GeneralRelativity.Constants.c "Speed of light";
  output Real deflection(unit = "rad") "Total deflection angle of a ray from infinity to infinity";
algorithm
  deflection := (1.0 + gamma) / 2.0 * 4.0 * GM / (c^2 * b);
  annotation(Documentation(info="<html>
<p>Leading-order deflection in the PPN metric. Half of the general-relativity value
<code>4GM/(c&sup2;b)</code> comes from <code>g<sub>00</sub></code> (the Newtonian
&quot;corpuscular&quot; value) and half from the spatial curvature <code>&gamma;</code>. A ray
grazing the Sun is bent by 1.75&Prime;.</p>
</html>"));
end ppnLightDeflection;
