within GeneralRelativity.Tensors;
function einstein "Einstein tensor G_ab = R_ab - 1/2 g_ab R"
  input Real g[4, 4] "Covariant metric";
  input Real Ric[4, 4] "Ricci tensor";
  input Real R "Ricci scalar";
  output Real G[4, 4] "G_ab";
algorithm
  G := Ric - 0.5 * R * g;
end einstein;
