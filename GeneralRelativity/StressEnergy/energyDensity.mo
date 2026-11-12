within GeneralRelativity.StressEnergy;
function energyDensity "Mass-energy density measured by an observer with four-velocity u: T_ab u^a u^b / c^2"
  input Real T[4, 4] "T_ab";
  input Real u[4] "Observer four-velocity, g(u,u) = -c^2";
  input Real c(unit = "m/s") = 1.0 "Speed of light";
  output Real rho(unit = "kg/m3");
algorithm
  rho := (u * (T * u)) / c^2;
  annotation(Inline = true);
end energyDensity;
