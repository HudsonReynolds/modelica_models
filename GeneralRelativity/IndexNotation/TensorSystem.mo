within GeneralRelativity.IndexNotation;
record TensorSystem "Metric and a table of named tensors at one event, the context an index expression is evaluated in"
  Real g[4, 4] "Covariant metric used to lower indices";
  Real gInv[4, 4] "Contravariant metric used to raise indices";
  Integer count "Number of table entries in use";
  String name[Internal.NameCapacity] "Tensor names";
  Integer rank[Internal.NameCapacity] "Number of indices of each entry";
  Integer offset[Internal.NameCapacity] "Index in data of each entry's first component";
  Integer used "Number of data entries in use";
  Real data[Internal.DataCapacity] "Fully covariant components of every entry, row-major with the last index fastest";
  annotation(Documentation(info="<html>
<p>Build one with
<a href=\"modelica://GeneralRelativity.IndexNotation.systemAt\">systemAt</a>,
<a href=\"modelica://GeneralRelativity.IndexNotation.systemFromCurvature\">systemFromCurvature</a>
or <a href=\"modelica://GeneralRelativity.IndexNotation.systemFromMetric\">systemFromMetric</a>
and extend it with the <code>define*</code> functions. Every tensor is stored with
all indices lowered; the evaluator raises what an expression asks for with
<code>gInv</code>. Two entries may share an offset, which is how the aliases
<code>Ric</code> and <code>Ricci</code> point at the data of <code>R</code> with
two indices.</p>
</html>"));
end TensorSystem;
