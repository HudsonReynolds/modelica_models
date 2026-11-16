within GeneralRelativity.IndexNotation.Internal;
function fixSlot "Take the component with slot k equal to a fixed value 1..4, dropping that slot"
  input GeneralRelativity.IndexNotation.TensorValue v;
  input Integer k "Slot, 1..v.rank";
  input Integer value "Array index 1..4 of the component to keep";
  output GeneralRelativity.IndexNotation.TensorValue r;
protected
  Integer nOut;
  Integer j;
  Integer o[MaxRank];
  Integer full[MaxRank];
algorithm
  r := emptyValue(v.rank - 1);
  j := 0;
  for q in 1:v.rank loop
    if q <> k then
      j := j + 1;
      r.index[j] := v.index[q];
      r.up[j] := v.up[q];
    end if;
  end for;
  nOut := power4(r.rank);
  for f in 1:nOut loop
    o := decode(f, r.rank);
    j := 0;
    for q in 1:v.rank loop
      if q == k then
        full[q] := value;
      else
        j := j + 1;
        full[q] := o[j];
      end if;
    end for;
    r.comp[f] := v.comp[encode(full, v.rank)];
  end for;
end fixSlot;
