within GeneralRelativity.Tensors;
function ricciScalar "Ricci scalar R = g^bd R_bd"
  input Real gInv[4, 4] "Contravariant metric";
  input Real Ric[4, 4] "Ricci tensor R_bd";
  output Real R "Ricci scalar";
algorithm
  R := sum(gInv[b, d] * Ric[b, d] for b in 1:4, d in 1:4);
end ricciScalar;
