within GeneralRelativity.IndexNotation.Internal;
function scale "Multiply every component of a value by a number"
  input GeneralRelativity.IndexNotation.TensorValue v;
  input Real factor;
  output GeneralRelativity.IndexNotation.TensorValue r;
algorithm
  r := v;
  r.comp := factor * v.comp;
end scale;
