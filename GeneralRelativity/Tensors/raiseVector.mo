within GeneralRelativity.Tensors;
function raiseVector "v^a = g^ab v_b"
  input Real gInv[4, 4] "Contravariant metric";
  input Real vLower[4] "Covariant components";
  output Real v[4] "Contravariant components";
algorithm
  v := gInv * vLower;
  annotation(Inline = true);
end raiseVector;
