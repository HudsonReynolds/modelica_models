within GeneralRelativity.Symmetry;
function rotationX "Rotation about the x axis in a spherical chart: -sin(phi) d_theta - cot(theta) cos(phi) d_phi"
  extends GeneralRelativity.Interfaces.VectorFieldFunction;
algorithm
  xi := {0.0, 0.0, -sin(x[4]), -cos(x[3]) / sin(x[3]) * cos(x[4])};
  annotation(Documentation(info="<html>
<p>The Cartesian field <code>y&part;<sub>z</sub> - z&part;<sub>y</sub></code> in the chart
<code>{ct, r, &theta;, &phi;}</code>. Together with
<a href=\"modelica://GeneralRelativity.Symmetry.rotationY\">rotationY</a> and
<a href=\"modelica://GeneralRelativity.Symmetry.axialRotation\">axialRotation</a>
it spans the rotation algebra: Killing for spherically symmetric metrics such
as Schwarzschild, not for Kerr. Singular on the polar axis. Parameters
<code>p</code> are ignored.</p>
</html>"));
end rotationX;
