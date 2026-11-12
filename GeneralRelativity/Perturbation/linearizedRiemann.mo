within GeneralRelativity.Perturbation;
function linearizedRiemann "Linear-order Riemann tensor from second derivatives of the perturbation"
  input Real ddh[4, 4, 4, 4] "ddh[k, l, a, b] = d^2 h_ab / (d x^k d x^l)";
  output Real R[4, 4, 4, 4] "R_abcd to first order in h";
algorithm
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        for d in 1:4 loop
          R[a, b, c, d] := 0.5 * (ddh[b, c, a, d] + ddh[a, d, b, c] - ddh[a, c, b, d] - ddh[b, d, a, c]);
        end for;
      end for;
    end for;
  end for;
end linearizedRiemann;
