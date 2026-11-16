within GeneralRelativity.IndexNotation;
function defineTensor3 "Add or replace a named rank-3 tensor; covariant components by default"
  input TensorSystem system;
  input String name;
  input Real A[4, 4, 4] "Components A[a, b, c]";
  input Boolean up[3] = {false, false, false} "Position of each index of A as passed, for example {true, false, false} for Gamma^a_bc";
  output TensorSystem result;
algorithm
  result := define(system, name, 3, Internal.flatten3(A), up);
end defineTensor3;
