within GeneralRelativity.StressEnergy;
function trace "T = g^ab T_ab; -rho c^2 + 3 p for a perfect fluid"
  input Real gInv[4, 4] "Contravariant metric";
  input Real T[4, 4] "T_ab";
  output Real result;
algorithm
  result := sum(gInv[a, b] * T[a, b] for a in 1:4, b in 1:4);
end trace;
