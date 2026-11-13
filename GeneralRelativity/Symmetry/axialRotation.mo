within GeneralRelativity.Symmetry;
function axialRotation "Rotation about the polar axis, d_phi, in a spherical chart {ct, r, theta, phi}"
  extends GeneralRelativity.Interfaces.VectorFieldFunction;
algorithm
  xi := {0.0, 0.0, 0.0, 1.0};
  annotation(Documentation(info="<html>
<p>The Cartesian rotation <code>x&part;<sub>y</sub> - y&part;<sub>x</sub></code>; Killing for
every axisymmetric metric in this chart (Schwarzschild, Kerr). Parameters
<code>p</code> are ignored.</p>
</html>"));
end axialRotation;
