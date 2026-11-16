within GeneralRelativity.IndexNotation;
record TensorValue "Value of an index expression: free index labels, their positions and flattened components"
  Integer rank "Number of free indices, at most Internal.MaxRank";
  String index[Internal.MaxRank] "Free index labels in slot order; unused slots hold the empty string";
  Boolean up[Internal.MaxRank] "true where the slot is contravariant (raised)";
  Real comp[Internal.Capacity] "Components in row-major order over the free slots, last slot fastest: slot values i_1..i_r are at comp[1 + sum_k (i_k - 1) 4^(r - k)]; entries beyond 4^rank are unused";
  annotation(Documentation(info="<html>
<p>Returned by <a href=\"modelica://GeneralRelativity.IndexNotation.evaluate\">evaluate</a>.
Use <a href=\"modelica://GeneralRelativity.IndexNotation.scalar\">scalar</a>,
<a href=\"modelica://GeneralRelativity.IndexNotation.vector\">vector</a>,
<a href=\"modelica://GeneralRelativity.IndexNotation.matrix\">matrix</a>,
<a href=\"modelica://GeneralRelativity.IndexNotation.tensor3\">tensor3</a> or
<a href=\"modelica://GeneralRelativity.IndexNotation.tensor4\">tensor4</a> to
unpack the components into an ordinary array, or
<a href=\"modelica://GeneralRelativity.IndexNotation.toString\">toString</a> to print them.</p>
</html>"));
end TensorValue;
