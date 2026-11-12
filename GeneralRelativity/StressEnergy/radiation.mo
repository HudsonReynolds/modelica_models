within GeneralRelativity.StressEnergy;
function radiation "Radiation fluid: perfect fluid with p = rho c^2 / 3, traceless"
  input Real g[4, 4] "Covariant metric";
  input Real u[4] "Contravariant four-velocity, g(u,u) = -c^2";
  input Real rho(unit = "kg/m3") "Energy density divided by c^2";
  input Real c(unit = "m/s") = 1.0 "Speed of light; 1 in geometric units";
  output Real T[4, 4] "T_ab";
algorithm
  T := perfectFluid(g, u, rho, rho * c^2 / 3.0, c);
end radiation;
