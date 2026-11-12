within GeneralRelativity.Tensors;
function innerProduct "g(u, v) = g_ab u^a v^b; negative for timelike pairs in the (-+++) signature"
  input Real g[4, 4] "Covariant metric";
  input Real u[4] "Contravariant components";
  input Real v[4] "Contravariant components";
  output Real result;
algorithm
  result := u * (g * v);
  annotation(Inline = true);
end innerProduct;
