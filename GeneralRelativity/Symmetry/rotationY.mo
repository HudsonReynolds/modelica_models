within GeneralRelativity.Symmetry;
function rotationY "Rotation about the y axis in a spherical chart: cos(phi) d_theta - cot(theta) sin(phi) d_phi"
  extends GeneralRelativity.Interfaces.VectorFieldFunction;
algorithm
  xi := {0.0, 0.0, cos(x[4]), -cos(x[3]) / sin(x[3]) * sin(x[4])};
  annotation(Documentation(info="<html>
<p>The Cartesian field <code>z&part;<sub>x</sub> - x&part;<sub>z</sub></code> in the chart
<code>{ct, r, &theta;, &phi;}</code>; see
<a href=\"modelica://GeneralRelativity.Symmetry.rotationX\">rotationX</a>.
Parameters <code>p</code> are ignored.</p>
</html>"));
end rotationY;
