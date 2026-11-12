within GeneralRelativity.Geodesics;
function acceleration "Geodesic acceleration a^i = -Gamma^i_jk u^j u^k"
  input Real Gamma[4, 4, 4] "Gamma[i, j, k] = Gamma^i_jk";
  input Real u[4] "Contravariant four-velocity dx/d lambda";
  output Real a[4] "d^2 x / d lambda^2";
algorithm
  for i in 1:4 loop
    a[i] := -sum(Gamma[i, j, k] * u[j] * u[k] for j in 1:4, k in 1:4);
  end for;
end acceleration;
