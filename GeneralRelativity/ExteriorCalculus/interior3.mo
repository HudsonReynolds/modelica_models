within GeneralRelativity.ExteriorCalculus;
function interior3 "Interior product of a vector with a 3-form, (i_X H)_bc = X^a H_abc"
  input Real X[4] "Vector X^a";
  input Real H[4, 4, 4] "3-form";
  output Real F[4, 4] "2-form";
algorithm
  for b in 1:4 loop
    for c in 1:4 loop
      F[b, c] := sum(X[a] * H[a, b, c] for a in 1:4);
    end for;
  end for;
end interior3;
