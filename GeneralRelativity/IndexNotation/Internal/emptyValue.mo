within GeneralRelativity.IndexNotation.Internal;
function emptyValue "TensorValue of the given rank with empty labels, lower positions and zero components"
  input Integer rank;
  output GeneralRelativity.IndexNotation.TensorValue v;
algorithm
  v.rank := rank;
  v.index := fill("", MaxRank);
  v.up := fill(false, MaxRank);
  v.comp := zeros(Capacity);
end emptyValue;
