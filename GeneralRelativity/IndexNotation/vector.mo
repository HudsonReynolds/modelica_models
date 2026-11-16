within GeneralRelativity.IndexNotation;
function vector "Evaluate an expression with one free index, such as \"R_{ab} u^b\", in the position the expression gives it"
  input TensorSystem system;
  input String expression;
  output Real value[4];
protected
  TensorValue v = evaluate(system, expression);
algorithm
  assert(v.rank == 1, "IndexNotation.vector: \"" + expression + "\" has free indices [" + Internal.indexList(v) + "], not one");
  value := v.comp[1:4];
end vector;
