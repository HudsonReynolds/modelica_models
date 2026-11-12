within GeneralRelativity.Metrics;
function cislunar "Weak-field Earth-Moon metric, Earth-centred inertial {ct, x, y, z} with the Moon on a circular orbit; p = {GM_earth, GM_moon, d, Omega, c, phi0}"
  extends GeneralRelativity.Interfaces.MetricFunction;
protected
  Real c = p[5];
  Real t = x[1] / c;
  Real moon[3] = GeneralRelativity.Perturbation.Potentials.circularBodyPosition(t, p[3], p[4], p[6]);
  Real phi = GeneralRelativity.Perturbation.Potentials.twoBody(x[2:4], {p[1], p[2], moon[1], moon[2], moon[3]});
algorithm
  g := GeneralRelativity.Perturbation.weakFieldMetric(phi, zeros(3), c);
  annotation(Documentation(info="<html>
<p>Superposed Newtonian potentials of Earth at the origin and the Moon at
<code>d (cos(&Omega;t + &phi;<sub>0</sub>), sin(&Omega;t + &phi;<sub>0</sub>), 0)</code> in the
weak-field metric. The metric is time dependent through <code>x[1] = ct</code>, so
the time derivative in <code>derivative</code> is nonzero but small
(<code>&Omega;/c</code> per metre). Equations of motion with first-order relativistic
corrections are in
<a href=\"modelica://GeneralRelativity.Examples.CislunarOrbit\">Examples.CislunarOrbit</a>.
The frame is inertial; the barycentre offset and lunar eccentricity are
neglected.</p>
</html>"));
end cislunar;
