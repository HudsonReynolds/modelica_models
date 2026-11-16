within GeneralRelativity.IndexNotation.Internal;
function decode "Multi-index (values 1..4) of a row-major flat index, last slot fastest; unused slots are 1"
  input Integer flat "1-based flat index";
  input Integer rank "Number of slots to decode";
  output Integer idx[MaxRank];
protected
  Integer f;
algorithm
  idx := fill(1, MaxRank);
  f := flat - 1;
  for k in rank:-1:1 loop
    idx[k] := mod(f, Dimension) + 1;
    f := div(f, Dimension);
  end for;
end decode;
