within GeneralRelativity.Geodesics;
function normResidual "g(u,u) + c^2 for a massive particle, g(u,u) for light; zero along an exact geodesic"
  input Real g[4, 4] "Covariant metric";
  input Real u[4] "Four-velocity";
  input Real c(unit = "m/s") = 1.0 "Speed of light; 0 for a null geodesic";
  output Real residual;
algorithm
  residual := u * (g * u) + c^2;
  annotation(Inline = true);
end normResidual;
