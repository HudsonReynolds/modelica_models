within GeneralRelativity.ExteriorCalculus;
function interior1 "Interior product of a vector with a 1-form, i_X w = X^a w_a"
  input Real X[4] "Vector X^a";
  input Real w[4] "1-form";
  output Real f "0-form";
algorithm
  f := X * w;
end interior1;
