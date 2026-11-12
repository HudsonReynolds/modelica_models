within GeneralRelativity.Tensors;
function christoffel "Christoffel symbols of the second kind Gamma[i, j, k] = Gamma^i_jk"
  input Real gInv[4, 4] "Contravariant metric";
  input Real dg[4, 4, 4] "dg[k, i, j] = d g_ij / d x^k";
  output Real Gamma[4, 4, 4] "Gamma[i, j, k] = Gamma^i_jk, symmetric in j, k";
protected
  Real lower[4, 4, 4] "lower[l, j, k] = Gamma_ljk";
algorithm
  for l in 1:4 loop
    for j in 1:4 loop
      for k in 1:4 loop
        lower[l, j, k] := 0.5 * (dg[j, l, k] + dg[k, l, j] - dg[l, j, k]);
      end for;
    end for;
  end for;
  for i in 1:4 loop
    for j in 1:4 loop
      for k in 1:4 loop
        Gamma[i, j, k] := sum(gInv[i, l] * lower[l, j, k] for l in 1:4);
      end for;
    end for;
  end for;
end christoffel;
