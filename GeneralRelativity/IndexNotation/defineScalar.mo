within GeneralRelativity.IndexNotation;
function defineScalar "Add or replace a named scalar"
  input TensorSystem system;
  input String name;
  input Real value;
  output TensorSystem result;
algorithm
  result := define(system, name, 0, {value});
end defineScalar;
