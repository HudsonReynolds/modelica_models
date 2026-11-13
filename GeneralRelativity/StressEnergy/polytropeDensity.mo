within GeneralRelativity.StressEnergy;
function polytropeDensity "Rest-mass density of a polytrope at given pressure, rho = (p/K)^(1/Gamma)"
  input Real pressure "Pressure (>= 0)";
  input Real K "Polytropic constant";
  input Real Gamma "Adiabatic index";
  output Real rho "Rest-mass density";
algorithm
  rho := (max(pressure, 0.0) / K)^(1.0 / Gamma);
  annotation(Documentation(info="<html>
<p>Inverse of
<a href=\"modelica://GeneralRelativity.StressEnergy.polytropePressure\">polytropePressure</a>,
used to close the TOV equations, which are integrated in the pressure. Negative
pressure (overshoot at the stellar surface) returns zero density.</p>
</html>"));
end polytropeDensity;
