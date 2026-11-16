within GeneralRelativity.IndexNotation;
function scalar "Evaluate an expression that has no free indices, such as \"R_{abcd} R^{abcd}\""
  input TensorSystem system;
  input String expression;
  output Real value;
protected
  TensorValue v = evaluate(system, expression);
algorithm
  assert(v.rank == 0, "IndexNotation.scalar: \"" + expression + "\" has free indices [" + Internal.indexList(v) + "]");
  value := v.comp[1];
end scalar;
