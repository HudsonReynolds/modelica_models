within GeneralRelativity.IndexNotation.Internal;
function register "Store fully covariant components under a name and rank, replacing an existing entry of the same name and rank"
  input GeneralRelativity.IndexNotation.TensorSystem sys;
  input String name;
  input Integer rank;
  input Real comp[:] "At least 4^rank components, row-major, last index fastest, all indices lowered";
  output GeneralRelativity.IndexNotation.TensorSystem result;
protected
  Integer n = power4(rank);
  Integer entry;
  Integer offset;
algorithm
  result := sys;
  assert(size(comp, 1) >= n, "IndexNotation: " + name + " with " + String(rank) + " indices needs " + String(n) + " components, got " + String(size(comp, 1)));
  entry := lookup(sys, name, rank);
  if entry == 0 then
    assert(sys.count < NameCapacity, "IndexNotation: the system already holds " + String(NameCapacity) + " names");
    assert(sys.used + n <= DataCapacity, "IndexNotation: no room for " + String(n) + " more components in the system");
    result.count := sys.count + 1;
    result.name[result.count] := name;
    result.rank[result.count] := rank;
    result.offset[result.count] := sys.used + 1;
    result.used := sys.used + n;
    entry := result.count;
  end if;
  offset := result.offset[entry];
  for i in 1:n loop
    result.data[offset + i - 1] := comp[i];
  end for;
end register;
