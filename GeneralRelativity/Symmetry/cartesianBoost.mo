within GeneralRelativity.Symmetry;
function cartesianBoost "Lorentz boost generator along the unit direction p = {px, py, pz} in a Cartesian chart {ct, x, y, z}"
  extends GeneralRelativity.Interfaces.VectorFieldFunction;
algorithm
  xi := {p[1] * x[2] + p[2] * x[3] + p[3] * x[4],
         p[1] * x[1],
         p[2] * x[1],
         p[3] * x[1]};
  annotation(Documentation(info="<html>
<p><code>&xi; = (p&middot;r) &part;<sub>0</sub> + x<sup>0</sup> p<sup>i</sup>&part;<sub>i</sub></code>; for
<code>p = {1, 0, 0}</code> this is <code>x&part;<sub>0</sub> + x<sup>0</sup>&part;<sub>x</sub></code>,
the generator of boosts along x. Its lowered form
<code>&xi;<sub>a</sub> = (-p&middot;r, x<sup>0</sup>p)</code> has an antisymmetric
gradient, so it is Killing in Minkowski space.</p>
</html>"));
end cartesianBoost;
