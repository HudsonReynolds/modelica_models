within GeneralRelativity.NewtonCartan;
function permutation "Three-dimensional Levi-Civita symbol epsilon_ijk with epsilon_123 = 1"
  input Integer i;
  input Integer j;
  input Integer k;
  output Real e;
algorithm
  e := (i - j) * (j - k) * (k - i) / 2.0;
  annotation(Inline = true);
end permutation;
