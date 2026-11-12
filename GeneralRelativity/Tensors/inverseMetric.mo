within GeneralRelativity.Tensors;
function inverseMetric "Contravariant metric g^ij by a pivoted dense solve"
  input Real g[4, 4] "Covariant metric";
  output Real gInv[4, 4] "Contravariant metric";
  output Boolean accepted "false when the metric is numerically singular";
algorithm
  (gInv, accepted) := LinearAlgebra.solve(g, identity(4));
  gInv := 0.5 * (gInv + transpose(gInv));
  annotation(Documentation(info="<html>
<p>The result is symmetrized. A false <code>accepted</code> flag means a
coordinate singularity such as the Schwarzschild horizon or the polar axis of a
spherical chart; curvature results built from it are not meaningful.</p>
</html>"));
end inverseMetric;
