within GeneralRelativity.IndexNotation.Internal;
function contractSlots "Sum over two slots of one value; when both have the same position the first is moved with the metric first"
  input GeneralRelativity.IndexNotation.TensorValue v;
  input Integer k "First slot";
  input Integer l "Second slot, different from k";
  input Real g[4, 4] "Covariant metric";
  input Real gInv[4, 4] "Contravariant metric";
  output GeneralRelativity.IndexNotation.TensorValue r;
protected
  GeneralRelativity.IndexNotation.TensorValue w;
  Integer nOut;
  Integer j;
  Integer o[MaxRank];
  Integer full[MaxRank];
  Real acc;
algorithm
  w := if v.up[k] == v.up[l] then setPosition(v, k, not v.up[k], g, gInv) else v;
  r := emptyValue(v.rank - 2);
  j := 0;
  for q in 1:v.rank loop
    if q <> k and q <> l then
      j := j + 1;
      r.index[j] := v.index[q];
      r.up[j] := w.up[q];
    end if;
  end for;
  nOut := power4(r.rank);
  for f in 1:nOut loop
    o := decode(f, r.rank);
    acc := 0.0;
    for m in 1:Dimension loop
      j := 0;
      for q in 1:v.rank loop
        if q == k or q == l then
          full[q] := m;
        else
          j := j + 1;
          full[q] := o[j];
        end if;
      end for;
      acc := acc + w.comp[encode(full, v.rank)];
    end for;
    r.comp[f] := acc;
  end for;
end contractSlots;
