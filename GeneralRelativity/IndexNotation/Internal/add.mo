within GeneralRelativity.IndexNotation.Internal;
function add "a + sign * b for values with the same free labels; b is permuted and raised or lowered to the slots and positions of a"
  input GeneralRelativity.IndexNotation.TensorValue a;
  input GeneralRelativity.IndexNotation.TensorValue b;
  input Real sign "+1 or -1";
  input Real g[4, 4] "Covariant metric";
  input Real gInv[4, 4] "Contravariant metric";
  output GeneralRelativity.IndexNotation.TensorValue r;
protected
  GeneralRelativity.IndexNotation.TensorValue bb;
  Integer perm[MaxRank] "Slot of b carrying the label of each slot of a";
  Integer n;
  Integer o[MaxRank];
  Integer ib[MaxRank];
algorithm
  assert(a.rank == b.rank,
    "IndexNotation: cannot add [" + indexList(a) + "] and [" + indexList(b) + "]: different numbers of free indices");
  bb := b;
  perm := zeros(MaxRank);
  for i in 1:a.rank loop
    for j in 1:b.rank loop
      if a.index[i] == b.index[j] then
        perm[i] := j;
      end if;
    end for;
    assert(perm[i] > 0,
      "IndexNotation: cannot add [" + indexList(a) + "] and [" + indexList(b) + "]: index " + a.index[i] + " is missing from the second term");
    if bb.up[perm[i]] <> a.up[i] then
      bb := setPosition(bb, perm[i], a.up[i], g, gInv);
    end if;
  end for;
  r := a;
  n := power4(a.rank);
  ib := fill(1, MaxRank);
  for f in 1:n loop
    o := decode(f, a.rank);
    for i in 1:a.rank loop
      ib[perm[i]] := o[i];
    end for;
    r.comp[f] := a.comp[f] + sign * bb.comp[encode(ib, a.rank)];
  end for;
end add;
