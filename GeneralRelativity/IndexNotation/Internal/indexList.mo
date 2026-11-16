within GeneralRelativity.IndexNotation.Internal;
function indexList "Free indices of a value as text, for example \"a^ b_\"; \"(scalar)\" for rank 0"
  input GeneralRelativity.IndexNotation.TensorValue v;
  output String text;
algorithm
  text := "";
  for k in 1:v.rank loop
    text := text + (if k > 1 then " " else "") + v.index[k] + (if v.up[k] then "^" else "_");
  end for;
  if v.rank == 0 then
    text := "(scalar)";
  end if;
end indexList;
