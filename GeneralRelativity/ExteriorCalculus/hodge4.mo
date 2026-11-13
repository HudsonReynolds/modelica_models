within GeneralRelativity.ExteriorCalculus;
function hodge4 "Hodge dual of a 4-form, *V = (1/4!) V^abcd eps_abcd"
  input Real V[4, 4, 4, 4] "4-form V_abcd (totally antisymmetric)";
  input Real g[4, 4] "Covariant metric";
  output Real f "0-form *V; *eps = -1 in Lorentzian signature";
protected
  Real gInv[4, 4];
  Boolean accepted;
  Real eps[4, 4, 4, 4];
  Real A[4, 4, 4, 4];
  Real B[4, 4, 4, 4];
algorithm
  (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(g);
  eps := leviCivitaTensor(g);
  // Raise one index at a time, alternating between two buffers.
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        for d in 1:4 loop
          A[a, b, c, d] := sum(gInv[a, m] * V[m, b, c, d] for m in 1:4);
        end for;
      end for;
    end for;
  end for;
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        for d in 1:4 loop
          B[a, b, c, d] := sum(gInv[b, m] * A[a, m, c, d] for m in 1:4);
        end for;
      end for;
    end for;
  end for;
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        for d in 1:4 loop
          A[a, b, c, d] := sum(gInv[c, m] * B[a, b, m, d] for m in 1:4);
        end for;
      end for;
    end for;
  end for;
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        for d in 1:4 loop
          B[a, b, c, d] := sum(gInv[d, m] * A[a, b, c, m] for m in 1:4);
        end for;
      end for;
    end for;
  end for;
  f := 0.0;
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        for d in 1:4 loop
          f := f + B[a, b, c, d] * eps[a, b, c, d] / 24.0;
        end for;
      end for;
    end for;
  end for;
end hodge4;
