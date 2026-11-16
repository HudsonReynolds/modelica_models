within GeneralRelativity.IndexNotation.Internal;
function encode "Row-major flat index (last slot fastest) of a multi-index with values 1..4"
  input Integer idx[:] "Slot values; only the first rank entries are read";
  input Integer rank "Number of slots";
  output Integer flat "1-based flat index";
algorithm
  flat := 0;
  for k in 1:rank loop
    flat := Dimension * flat + idx[k] - 1;
  end for;
  flat := flat + 1;
end encode;
