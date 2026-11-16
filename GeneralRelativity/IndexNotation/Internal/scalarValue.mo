within GeneralRelativity.IndexNotation.Internal;
function scalarValue "Rank-0 TensorValue holding one number"
  input Real value;
  output GeneralRelativity.IndexNotation.TensorValue v;
algorithm
  v := emptyValue(0);
  v.comp[1] := value;
end scalarValue;
