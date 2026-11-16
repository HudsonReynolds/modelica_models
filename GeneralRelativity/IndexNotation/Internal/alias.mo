within GeneralRelativity.IndexNotation.Internal;
function alias "Add a second name for an existing entry; both names then read the same components"
  input GeneralRelativity.IndexNotation.TensorSystem sys;
  input String name "New name";
  input String target "Existing name";
  input Integer rank "Number of indices of both";
  output GeneralRelativity.IndexNotation.TensorSystem result;
protected
  Integer source = lookup(sys, target, rank);
  Integer entry = lookup(sys, name, rank);
algorithm
  result := sys;
  assert(source > 0, "IndexNotation: cannot alias " + name + " to the undefined tensor " + target + " with " + String(rank) + " indices");
  if entry == 0 then
    assert(sys.count < NameCapacity, "IndexNotation: the system already holds " + String(NameCapacity) + " names");
    result.count := sys.count + 1;
    entry := result.count;
    result.name[entry] := name;
    result.rank[entry] := rank;
  end if;
  result.offset[entry] := sys.offset[source];
end alias;
