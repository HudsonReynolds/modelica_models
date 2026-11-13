within GeneralRelativity.ExteriorCalculus;
function hodge3 "Hodge dual of a 3-form, (*H)_d = (1/3!) H^abc eps_abcd"
  input Real H[4, 4, 4] "3-form H_abc (totally antisymmetric)";
  input Real g[4, 4] "Covariant metric";
  output Real w[4] "1-form *H";
protected
  Real gInv[4, 4];
  Boolean accepted;
  Real eps[4, 4, 4, 4];
  Real T1[4, 4, 4];
  Real T2[4, 4, 4];
  Real HUp[4, 4, 4];
algorithm
  (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(g);
  eps := leviCivitaTensor(g);
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        T1[a, b, c] := sum(gInv[a, d] * H[d, b, c] for d in 1:4);
      end for;
    end for;
  end for;
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        T2[a, b, c] := sum(gInv[b, d] * T1[a, d, c] for d in 1:4);
      end for;
    end for;
  end for;
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        HUp[a, b, c] := sum(gInv[c, d] * T2[a, b, d] for d in 1:4);
      end for;
    end for;
  end for;
  for d in 1:4 loop
    w[d] := 0.0;
    for a in 1:4 loop
      for b in 1:4 loop
        for c in 1:4 loop
          w[d] := w[d] + HUp[a, b, c] * eps[a, b, c, d] / 6.0;
        end for;
      end for;
    end for;
  end for;
end hodge3;
