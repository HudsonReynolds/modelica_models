within GeneralRelativity.Symmetry;
function cartesianRotation "Spatial rotation p x r about the axis p = {px, py, pz} in a Cartesian chart {ct, x, y, z}"
  extends GeneralRelativity.Interfaces.VectorFieldFunction;
algorithm
  xi := {0.0,
         p[2] * x[4] - p[3] * x[3],
         p[3] * x[2] - p[1] * x[4],
         p[1] * x[3] - p[2] * x[2]};
  annotation(Documentation(info="<html>
<p><code>p = {0, 0, 1}</code> gives <code>x&part;<sub>y</sub> - y&part;<sub>x</sub></code>,
a counter-clockwise rotation about the z axis. Killing in Minkowski space.</p>
</html>"));
end cartesianRotation;
