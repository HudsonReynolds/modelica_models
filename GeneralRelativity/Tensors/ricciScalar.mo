within GeneralRelativity.Tensors;
function ricciScalar "Ricci scalar R = g^bd R_bd"
  input Real gInv[:, size(gInv, 1)] "Contravariant metric";
  input Real Ric[size(gInv, 1), size(gInv, 1)] "Ricci tensor R_bd";
  output Real R "Ricci scalar";
protected
  Integer n = size(gInv, 1);
algorithm
  R := sum(gInv[b, d] * Ric[b, d] for b in 1:n, d in 1:n);
end ricciScalar;
