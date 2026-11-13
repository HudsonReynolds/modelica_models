within GeneralRelativity.Tensors;
function inverseMetric "Contravariant metric g^ij by a pivoted dense solve, in any dimension"
  input Real g[:, size(g, 1)] "Covariant metric";
  output Real gInv[size(g, 1), size(g, 1)] "Contravariant metric";
  output Boolean accepted "false when the metric is numerically singular";
algorithm
  (gInv, accepted) := LinearAlgebra.solve(g, identity(size(g, 1)));
  gInv := 0.5 * (gInv + transpose(gInv));
  annotation(Documentation(info="<html>
<p>The result is symmetrized. A false <code>accepted</code> flag means a
coordinate singularity such as the Schwarzschild horizon or the polar axis of a
spherical chart; curvature results built from it are not meaningful.</p>
</html>"));
end inverseMetric;
