within GeneralRelativity.StressEnergy;
function polytropePressure "Polytropic equation of state p = K rho^Gamma"
  input Real rho "Rest-mass density (>= 0)";
  input Real K "Polytropic constant, in units with p and rho consistent";
  input Real Gamma "Adiabatic index Gamma = 1 + 1/n";
  output Real pressure "Pressure";
algorithm
  pressure := K * max(rho, 0.0)^Gamma;
  annotation(Documentation(info="<html>
<p>The polytropic index is <code>n = 1/(&Gamma; - 1)</code>; <code>&Gamma; = 2</code>
(<code>n = 1</code>) is the textbook neutron-star model and <code>&Gamma; = 5/3</code>
(<code>n = 3/2</code>) a non-relativistic degenerate gas. Negative densities are
clipped to zero so the function is safe at a stellar surface. Units: with
<code>c = G = 1</code>, <code>K</code> has the dimension length<sup>2/n</sup>. See
<a href=\"modelica://GeneralRelativity.StressEnergy.polytropeEnergyDensity\">polytropeEnergyDensity</a>
and <a href=\"modelica://GeneralRelativity.StressEnergy.polytropeDensity\">polytropeDensity</a>.</p>
</html>"));
end polytropePressure;
