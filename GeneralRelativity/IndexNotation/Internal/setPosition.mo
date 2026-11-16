within GeneralRelativity.IndexNotation.Internal;
function setPosition "Raise or lower slot k of a value with the metric so that it has the requested position; no-op when it already does"
  input GeneralRelativity.IndexNotation.TensorValue v;
  input Integer k "Slot, 1..v.rank";
  input Boolean toUp "true to raise, false to lower";
  input Real g[4, 4] "Covariant metric";
  input Real gInv[4, 4] "Contravariant metric";
  output GeneralRelativity.IndexNotation.TensorValue r;
protected
  Integer n = power4(v.rank);
  Integer s = power4(v.rank - k) "Stride of slot k";
  Integer i;
  Integer base;
  Real M[4, 4] = if toUp then gInv else g;
algorithm
  r := v;
  if v.up[k] <> toUp then
    r.up[k] := toUp;
    for f in 1:n loop
      i := mod(div(f - 1, s), Dimension) + 1;
      base := f - (i - 1) * s;
      r.comp[f] := sum(M[i, m] * v.comp[base + (m - 1) * s] for m in 1:Dimension);
    end for;
  end if;
end setPosition;
