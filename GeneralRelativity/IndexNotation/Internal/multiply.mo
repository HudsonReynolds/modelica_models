within GeneralRelativity.IndexNotation.Internal;
function multiply "Outer product of two values with every label they share summed; free slots of the first value come first"
  input GeneralRelativity.IndexNotation.TensorValue a;
  input GeneralRelativity.IndexNotation.TensorValue b;
  input Real g[4, 4] "Covariant metric";
  input Real gInv[4, 4] "Contravariant metric";
  output GeneralRelativity.IndexNotation.TensorValue r;
protected
  GeneralRelativity.IndexNotation.TensorValue aa;
  Integer m "Number of contracted pairs";
  Integer pairA[MaxRank] "Pair number of each slot of a, 0 when free";
  Integer pairB[MaxRank];
  Integer freeA[MaxRank] "Result slot of each free slot of a";
  Integer freeB[MaxRank];
  Integer partner[MaxRank] "Slot of b paired with each slot of a";
  Integer nf;
  Integer nOut;
  Integer nC;
  Integer o[MaxRank];
  Integer cidx[MaxRank];
  Integer ia[MaxRank];
  Integer ib[MaxRank];
  Real acc;
algorithm
  m := 0;
  pairA := zeros(MaxRank);
  pairB := zeros(MaxRank);
  partner := zeros(MaxRank);
  for i in 1:a.rank loop
    for j in 1:b.rank loop
      if a.index[i] == b.index[j] then
        m := m + 1;
        pairA[i] := m;
        pairB[j] := m;
        partner[i] := j;
      end if;
    end for;
  end for;
  assert(a.rank + b.rank - 2 * m <= MaxRank,
    "IndexNotation: the product of [" + indexList(a) + "] and [" + indexList(b) + "] has more than "
    + String(MaxRank) + " free indices");
  // A contracted pair needs one upper and one lower index: move the slot of a to the opposite of b's.
  aa := a;
  for i in 1:a.rank loop
    if pairA[i] > 0 and aa.up[i] == b.up[partner[i]] then
      aa := setPosition(aa, i, not b.up[partner[i]], g, gInv);
    end if;
  end for;
  r := emptyValue(a.rank + b.rank - 2 * m);
  nf := 0;
  freeA := zeros(MaxRank);
  freeB := zeros(MaxRank);
  for i in 1:a.rank loop
    if pairA[i] == 0 then
      nf := nf + 1;
      r.index[nf] := aa.index[i];
      r.up[nf] := aa.up[i];
      freeA[i] := nf;
    end if;
  end for;
  for j in 1:b.rank loop
    if pairB[j] == 0 then
      nf := nf + 1;
      r.index[nf] := b.index[j];
      r.up[nf] := b.up[j];
      freeB[j] := nf;
    end if;
  end for;
  nOut := power4(r.rank);
  nC := power4(m);
  ia := fill(1, MaxRank);
  ib := fill(1, MaxRank);
  for f in 1:nOut loop
    o := decode(f, r.rank);
    acc := 0.0;
    for h in 1:nC loop
      cidx := decode(h, m);
      for i in 1:a.rank loop
        ia[i] := if pairA[i] > 0 then cidx[pairA[i]] else o[freeA[i]];
      end for;
      for j in 1:b.rank loop
        ib[j] := if pairB[j] > 0 then cidx[pairB[j]] else o[freeB[j]];
      end for;
      acc := acc + aa.comp[encode(ia, a.rank)] * b.comp[encode(ib, b.rank)];
    end for;
    r.comp[f] := acc;
  end for;
end multiply;
