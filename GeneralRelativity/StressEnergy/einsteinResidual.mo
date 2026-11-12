within GeneralRelativity.StressEnergy;
function einsteinResidual "Field-equation residual G_ab + Lambda g_ab - 8 pi G / c^4 T_ab"
  input Real G[4, 4] "Einstein tensor";
  input Real g[4, 4] "Covariant metric";
  input Real T[4, 4] "Stress-energy T_ab";
  input Real gravitationalConstant(unit = "m3/(kg.s2)") = 1.0 "Newton's constant";
  input Real c(unit = "m/s") = 1.0 "Speed of light";
  input Real Lambda(unit = "1/m2") = 0.0 "Cosmological constant";
  output Real residual[4, 4] "Zero when the metric solves the field equations with this source";
algorithm
  residual := G + Lambda * g - 8.0 * GeneralRelativity.Constants.pi * gravitationalConstant / c^4 * T;
end einsteinResidual;
