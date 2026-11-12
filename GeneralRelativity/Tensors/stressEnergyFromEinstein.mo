within GeneralRelativity.Tensors;
function stressEnergyFromEinstein "Stress-energy tensor implied by a metric: T_ab = c^4/(8 pi G) (G_ab + Lambda g_ab)"
  input Real G[4, 4] "Einstein tensor G_ab";
  input Real g[4, 4] "Covariant metric";
  input Real gravitationalConstant(unit = "m3/(kg.s2)") = 1.0 "Newton's constant G; 1 in geometric units";
  input Real c(unit = "m/s") = 1.0 "Speed of light; 1 in geometric units";
  input Real Lambda(unit = "1/m2") = 0.0 "Cosmological constant";
  output Real T[4, 4] "T_ab; T[1, 1] is the energy density for a static observer in an orthonormal time direction";
algorithm
  T := c^4 / (8.0 * GeneralRelativity.Constants.pi * gravitationalConstant) * (G + Lambda * g);
  annotation(Documentation(info="<html>
<p>Reading the field equations backwards: a metric that is a vacuum solution
returns zero; the weak-field Earth metric returns the Newtonian mass density
through <code>T<sub>11</sub> = &rho;c&sup2;</code>. Compare with the source you
intended using
<a href=\"modelica://GeneralRelativity.StressEnergy.einsteinResidual\">StressEnergy.einsteinResidual</a>.</p>
</html>"));
end stressEnergyFromEinstein;
