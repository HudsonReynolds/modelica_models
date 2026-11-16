within GeneralRelativity.IndexNotation;
package Internal "Index algebra on flattened components and the recursive-descent parser behind evaluate; not part of the public API"
  constant Integer Dimension = 4 "Spacetime dimension; every index runs 1..4";
  constant Integer MaxRank = 6 "Largest number of free indices an intermediate value may carry";
  constant Integer Capacity = 4096 "Dimension^MaxRank, the component storage of a TensorValue";
  constant Integer NameCapacity = 64 "Table entries a TensorSystem can hold";
  constant Integer DataCapacity = 8192 "Components a TensorSystem can hold; the built-in curvature set uses 917";
  annotation(Documentation(info="<html>
<p>Values are <a href=\"modelica://GeneralRelativity.IndexNotation.TensorValue\">TensorValue</a>
records: the components of a rank-<code>r</code> tensor are stored row-major
(last index fastest) in the first <code>4^r</code> entries of a fixed array, see
<a href=\"modelica://GeneralRelativity.IndexNotation.Internal.encode\">encode</a> and
<a href=\"modelica://GeneralRelativity.IndexNotation.Internal.decode\">decode</a>.
The algebra is
<a href=\"modelica://GeneralRelativity.IndexNotation.Internal.setPosition\">setPosition</a>
(raise or lower one slot),
<a href=\"modelica://GeneralRelativity.IndexNotation.Internal.contractSlots\">contractSlots</a>,
<a href=\"modelica://GeneralRelativity.IndexNotation.Internal.fixSlot\">fixSlot</a>,
<a href=\"modelica://GeneralRelativity.IndexNotation.Internal.multiply\">multiply</a>
(outer product with contraction of shared labels),
<a href=\"modelica://GeneralRelativity.IndexNotation.Internal.add\">add</a> and
<a href=\"modelica://GeneralRelativity.IndexNotation.Internal.scale\">scale</a>.
The parser is the mutually recursive trio
<a href=\"modelica://GeneralRelativity.IndexNotation.Internal.parseExpression\">parseExpression</a>
(sums), <a href=\"modelica://GeneralRelativity.IndexNotation.Internal.parseTerm\">parseTerm</a>
(products and quotients) and
<a href=\"modelica://GeneralRelativity.IndexNotation.Internal.parseFactor\">parseFactor</a>
(numbers, parentheses and
<a href=\"modelica://GeneralRelativity.IndexNotation.Internal.parseTensor\">tensor references</a>).</p>
</html>"));
end Internal;
