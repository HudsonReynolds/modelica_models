within GeneralRelativity.Tensors;
function einstein "Einstein tensor G_ab = R_ab - 1/2 g_ab R"
  input Real g[:, size(g, 1)] "Covariant metric";
  input Real Ric[size(g, 1), size(g, 1)] "Ricci tensor";
  input Real R "Ricci scalar";
  output Real G[size(g, 1), size(g, 1)] "G_ab";
algorithm
  G := Ric - 0.5 * R * g;
end einstein;
