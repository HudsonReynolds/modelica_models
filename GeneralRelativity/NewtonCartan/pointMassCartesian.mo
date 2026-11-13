within GeneralRelativity.NewtonCartan;
function pointMassCartesian "Newton-Cartan connection of a point mass in a Cartesian chart, optionally rotating; p = {rs} or {rs, wx, wy, wz}"
  extends GeneralRelativity.Interfaces.ConnectionFunction;
protected
  Real rs = p[1] "2GM/c^2";
  Real w[3] = if size(p, 1) >= 4 then p[2:4] else zeros(3) "Chart angular velocity Omega/c";
  Real position[3] = x[2:4];
  Real r = sqrt(position * position);
  Real centrifugal[3];
algorithm
  Gamma := zeros(4, 4, 4);
  centrifugal := cross(w, cross(w, position));
  for i in 1:3 loop
    Gamma[i + 1, 1, 1] := (if rs > 0.0 then rs * position[i] / (2.0 * r^3) else 0.0) + centrifugal[i];
    for j in 1:3 loop
      // Coriolis: Gamma^i_0j v^j = (w x v)^i, antisymmetric in i, j.
      Gamma[i + 1, 1, j + 1] := sum(GeneralRelativity.NewtonCartan.permutation(i, k, j) * w[k] for k in 1:3);
      Gamma[i + 1, j + 1, 1] := Gamma[i + 1, 1, j + 1];
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>With <code>&lambda; = ct</code> and <code>w = &Omega;/c</code>, the geodesic equation reads
<code>a = &minus;GM x/r&sup3; &minus; 2&Omega; &times; v &minus; &Omega; &times; (&Omega; &times; x)</code>.
The Coriolis coefficients <code>&Gamma;<sup>i</sup><sub>0j</sub></code> are antisymmetric in
<code>i, j</code>, which is what keeps <code>&nabla;h = 0</code>. With <code>r<sub>s</sub> = 0</code>
the connection is flat however fast the chart rotates: inertial forces are
coordinate effects; tidal curvature is not.</p>
</html>"));
end pointMassCartesian;
