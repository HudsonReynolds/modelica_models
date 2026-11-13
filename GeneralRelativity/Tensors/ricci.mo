within GeneralRelativity.Tensors;
function ricci "Ricci tensor R_bd = g^ac R_abcd"
  input Real gInv[:, size(gInv, 1)] "Contravariant metric";
  input Real RLower[size(gInv, 1), size(gInv, 1), size(gInv, 1), size(gInv, 1)] "R_abcd";
  output Real Ric[size(gInv, 1), size(gInv, 1)] "R_bd";
protected
  Integer n = size(gInv, 1);
algorithm
  for b in 1:n loop
    for d in 1:n loop
      Ric[b, d] := sum(gInv[a, c] * RLower[a, b, c, d] for a in 1:n, c in 1:n);
    end for;
  end for;
end ricci;
