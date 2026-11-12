within GeneralRelativity.Tensors;
function lowerVector "v_a = g_ab v^b"
  input Real g[4, 4] "Covariant metric";
  input Real v[4] "Contravariant components";
  output Real vLower[4] "Covariant components";
algorithm
  vLower := g * v;
  annotation(Inline = true);
end lowerVector;
