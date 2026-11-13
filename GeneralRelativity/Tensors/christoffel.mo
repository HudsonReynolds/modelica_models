within GeneralRelativity.Tensors;
function christoffel "Christoffel symbols of the second kind Gamma[i, j, k] = Gamma^i_jk, in any dimension"
  input Real gInv[:, size(gInv, 1)] "Contravariant metric";
  input Real dg[size(gInv, 1), size(gInv, 1), size(gInv, 1)] "dg[k, i, j] = d g_ij / d x^k";
  output Real Gamma[size(gInv, 1), size(gInv, 1), size(gInv, 1)] "Gamma[i, j, k] = Gamma^i_jk, symmetric in j, k";
protected
  Integer n = size(gInv, 1);
  Real lower[size(gInv, 1), size(gInv, 1), size(gInv, 1)] "lower[l, j, k] = Gamma_ljk";
algorithm
  for l in 1:n loop
    for j in 1:n loop
      for k in 1:n loop
        lower[l, j, k] := 0.5 * (dg[j, l, k] + dg[k, l, j] - dg[l, j, k]);
      end for;
    end for;
  end for;
  for i in 1:n loop
    for j in 1:n loop
      for k in 1:n loop
        Gamma[i, j, k] := sum(gInv[i, l] * lower[l, j, k] for l in 1:n);
      end for;
    end for;
  end for;
end christoffel;
