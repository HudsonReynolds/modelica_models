within GeneralRelativity.Connections;
function torsion "Torsion T^i_jk = Gamma^i_jk - Gamma^i_kj"
  input Real Gamma[4, 4, 4] "Gamma[i, j, k] = Gamma^i_jk";
  output Real T[4, 4, 4] "T^i_jk, antisymmetric in j, k";
algorithm
  for i in 1:4 loop
    for j in 1:4 loop
      for k in 1:4 loop
        T[i, j, k] := Gamma[i, j, k] - Gamma[i, k, j];
      end for;
    end for;
  end for;
end torsion;
