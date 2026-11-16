within GeneralRelativity.IndexNotation;
function defineTensor4 "Add or replace a named rank-4 tensor; covariant components by default"
  input TensorSystem system;
  input String name;
  input Real A[4, 4, 4, 4] "Components A[a, b, c, d]";
  input Boolean up[4] = {false, false, false, false} "Position of each index of A as passed";
  output TensorSystem result;
algorithm
  result := define(system, name, 4, Internal.flatten4(A), up);
end defineTensor4;
