within GeneralRelativity.Perturbation;
function newtonianPotential "Potential read from g_00: Phi = -c^2 (1 + g_00) / 2"
  input Real g[4, 4] "Covariant metric in a nearly Cartesian chart";
  input Real c(unit = "m/s") = 1.0 "Speed of light";
  output Real phi(unit = "m2/s2");
algorithm
  phi := -0.5 * c^2 * (1.0 + g[1, 1]);
  annotation(Inline = true);
end newtonianPotential;
