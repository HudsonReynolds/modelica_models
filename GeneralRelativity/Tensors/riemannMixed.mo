within GeneralRelativity.Tensors;
function riemannMixed "Riemann tensor with the first index raised, R[a, b, c, d] = R^a_bcd"
  input Real gInv[4, 4] "Contravariant metric";
  input Real RLower[4, 4, 4, 4] "R_abcd";
  output Real R[4, 4, 4, 4] "R^a_bcd";
algorithm
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        for d in 1:4 loop
          R[a, b, c, d] := sum(gInv[a, e] * RLower[e, b, c, d] for e in 1:4);
        end for;
      end for;
    end for;
  end for;
end riemannMixed;
