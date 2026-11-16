within GeneralRelativity.IndexNotation;
function defineVector "Add or replace a named vector; contravariant components by default"
  input TensorSystem system;
  input String name;
  input Real v[4] "Components";
  input Boolean up = true "true when v holds contravariant components v^a, false for covariant v_a";
  output TensorSystem result;
algorithm
  result := define(system, name, 1, v, {up});
end defineVector;
