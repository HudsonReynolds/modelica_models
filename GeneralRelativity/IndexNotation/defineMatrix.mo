within GeneralRelativity.IndexNotation;
function defineMatrix "Add or replace a named rank-2 tensor; covariant components by default"
  input TensorSystem system;
  input String name;
  input Real M[4, 4] "Components M[a, b]";
  input Boolean up[2] = {false, false} "Position of each index of M as passed; {true, true} for M^ab, {true, false} for M^a_b";
  output TensorSystem result;
algorithm
  result := define(system, name, 2, Internal.flatten2(M), up);
end defineMatrix;
