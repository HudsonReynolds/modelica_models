within GeneralRelativity.Tensors;
function riemannMixed "Riemann tensor with the first index raised, R[a, b, c, d] = R^a_bcd"
  input Real gInv[:, size(gInv, 1)] "Contravariant metric";
  input Real RLower[size(gInv, 1), size(gInv, 1), size(gInv, 1), size(gInv, 1)] "R_abcd";
  output Real R[size(gInv, 1), size(gInv, 1), size(gInv, 1), size(gInv, 1)] "R^a_bcd";
protected
  Integer n = size(gInv, 1);
algorithm
  for a in 1:n loop
    for b in 1:n loop
      for c in 1:n loop
        for d in 1:n loop
          R[a, b, c, d] := sum(gInv[a, e] * RLower[e, b, c, d] for e in 1:n);
        end for;
      end for;
    end for;
  end for;
end riemannMixed;
