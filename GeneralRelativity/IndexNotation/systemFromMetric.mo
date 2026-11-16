within GeneralRelativity.IndexNotation;
function systemFromMetric "System holding only a metric: g, delta and the Levi-Civita tensor epsilon, for algebra on user-defined tensors"
  input Real g[4, 4] "Covariant metric";
  output TensorSystem sys;
protected
  Boolean accepted;
algorithm
  sys.g := g;
  (sys.gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(g);
  assert(accepted, "IndexNotation.systemFromMetric: the metric is numerically singular");
  sys.count := 0;
  sys.name := fill("", Internal.NameCapacity);
  sys.rank := zeros(Internal.NameCapacity);
  sys.offset := zeros(Internal.NameCapacity);
  sys.used := 0;
  sys.data := zeros(Internal.DataCapacity);
  sys := Internal.register(sys, "g", 2, Internal.flatten2(g));
  sys := Internal.alias(sys, "delta", "g", 2);
  sys := Internal.register(sys, "epsilon", 4, Internal.flatten4(GeneralRelativity.ExteriorCalculus.leviCivitaTensor(g)));
  sys := Internal.alias(sys, "eps", "epsilon", 4);
  annotation(Documentation(info="<html>
<p><code>delta</code> is the metric under another name: with one index up and one
down it is the Kronecker delta. <code>epsilon_{0123} = +sqrt(|det g|)</code>, see
<a href=\"modelica://GeneralRelativity.ExteriorCalculus.leviCivitaTensor\">ExteriorCalculus.leviCivitaTensor</a>.</p>
</html>"));
end systemFromMetric;
