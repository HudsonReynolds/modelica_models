within GeneralRelativity.IndexNotation.Internal;
function power4 "4^e as an Integer for e >= 0; the number of components of a rank-e tensor"
  input Integer e;
  output Integer p;
algorithm
  p := 1;
  for i in 1:e loop
    p := Dimension * p;
  end for;
  annotation(Inline = true);
end power4;
