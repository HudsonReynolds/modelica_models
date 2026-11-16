within GeneralRelativity.IndexNotation.Internal;
function lookup "Table entry of a name with a given number of indices, 0 when the system has none"
  input GeneralRelativity.IndexNotation.TensorSystem sys;
  input String name;
  input Integer rank;
  output Integer entry;
algorithm
  entry := 0;
  for i in 1:sys.count loop
    if entry == 0 and sys.rank[i] == rank and sys.name[i] == name then
      entry := i;
    end if;
  end for;
end lookup;
