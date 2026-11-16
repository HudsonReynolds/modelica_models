within GeneralRelativity.IndexNotation;
function tensor4 "Evaluate an expression with four free indices, such as \"R^{ab}_{cd}\", as A[first, second, third, fourth index]"
  input TensorSystem system;
  input String expression;
  output Real value[4, 4, 4, 4];
protected
  TensorValue v = evaluate(system, expression);
algorithm
  assert(v.rank == 4, "IndexNotation.tensor4: \"" + expression + "\" has free indices [" + Internal.indexList(v) + "], not four");
  for i in 1:4 loop
    for j in 1:4 loop
      for k in 1:4 loop
        for l in 1:4 loop
          value[i, j, k, l] := v.comp[(((i - 1) * 4 + (j - 1)) * 4 + (k - 1)) * 4 + l];
        end for;
      end for;
    end for;
  end for;
end tensor4;
