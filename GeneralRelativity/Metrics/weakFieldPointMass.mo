within GeneralRelativity.Metrics;
function weakFieldPointMass "Linearized point-mass metric in isotropic Cartesian coordinates {ct, x, y, z}; p = {GM, c}"
  extends GeneralRelativity.Interfaces.MetricFunction;
protected
  Real phi = GeneralRelativity.Perturbation.Potentials.pointMass(x[2:4], {p[1]});
algorithm
  g := GeneralRelativity.Perturbation.weakFieldMetric(phi, zeros(3), p[2]);
  annotation(Documentation(info="<html>
<p><code>g = &eta; + h</code> with <code>h<sub>00</sub> = h<sub>ii</sub> = -2&Phi;/c&sup2;</code>
and <code>&Phi; = -GM/r</code>. This is Schwarzschild in isotropic coordinates to first
order in <code>GM/(c&sup2;r)</code>.</p>
</html>"));
end weakFieldPointMass;
