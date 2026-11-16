within GeneralRelativity.IndexNotation;
function matrix "Evaluate an expression with two free indices, such as \"R^a_{bac}\", as M[first index, second index]"
  input TensorSystem system;
  input String expression;
  output Real value[4, 4];
protected
  TensorValue v = evaluate(system, expression);
algorithm
  assert(v.rank == 2, "IndexNotation.matrix: \"" + expression + "\" has free indices [" + Internal.indexList(v) + "], not two");
  for i in 1:4 loop
    for j in 1:4 loop
      value[i, j] := v.comp[(i - 1) * 4 + j];
    end for;
  end for;
end matrix;
