within GeneralRelativity.Perturbation;
function linearizedRicci "Linear-order Ricci tensor R_bd = eta^ac R_abcd"
  input Real R[4, 4, 4, 4] "Linearized R_abcd";
  output Real Ric[4, 4];
protected
  constant Real eta[4] = {-1.0, 1.0, 1.0, 1.0};
algorithm
  for b in 1:4 loop
    for d in 1:4 loop
      Ric[b, d] := sum(eta[a] * R[a, b, a, d] for a in 1:4);
    end for;
  end for;
end linearizedRicci;
