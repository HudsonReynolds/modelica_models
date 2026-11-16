within GeneralRelativity.IndexNotation;
function tensor3 "Evaluate an expression with three free indices, such as \"Gamma^i_{jk}\", as A[first, second, third index]"
  input TensorSystem system;
  input String expression;
  output Real value[4, 4, 4];
protected
  TensorValue v = evaluate(system, expression);
algorithm
  assert(v.rank == 3, "IndexNotation.tensor3: \"" + expression + "\" has free indices [" + Internal.indexList(v) + "], not three");
  for i in 1:4 loop
    for j in 1:4 loop
      for k in 1:4 loop
        value[i, j, k] := v.comp[((i - 1) * 4 + (j - 1)) * 4 + k];
      end for;
    end for;
  end for;
end tensor3;
