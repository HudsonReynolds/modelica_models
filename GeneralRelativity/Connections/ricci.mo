within GeneralRelativity.Connections;
function ricci "Ricci tensor R_sn = R^r_srn of a connection; symmetric only for a torsion-free metric connection"
  input Real R[4, 4, 4, 4] "R^r_smn";
  output Real Ric[4, 4] "R_sn";
algorithm
  for s in 1:4 loop
    for n in 1:4 loop
      Ric[s, n] := sum(R[r, s, r, n] for r in 1:4);
    end for;
  end for;
end ricci;
