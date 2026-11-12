within GeneralRelativity.Tensors;
function ricci "Ricci tensor R_bd = g^ac R_abcd"
  input Real gInv[4, 4] "Contravariant metric";
  input Real RLower[4, 4, 4, 4] "R_abcd";
  output Real Ric[4, 4] "R_bd";
algorithm
  for b in 1:4 loop
    for d in 1:4 loop
      Ric[b, d] := sum(gInv[a, c] * RLower[a, b, c, d] for a in 1:4, c in 1:4);
    end for;
  end for;
end ricci;
