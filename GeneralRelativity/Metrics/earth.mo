within GeneralRelativity.Metrics;
function earth "Weak-field Earth metric with J2 oblateness and Lense-Thirring frame dragging; p = {GM, R, J2, c, J}"
  extends GeneralRelativity.Interfaces.MetricFunction;
protected
  Real GM = p[1] "Gravitational parameter";
  Real R = p[2] "Equatorial radius";
  Real J2 = p[3] "Second zonal harmonic";
  Real c = p[4] "Speed of light";
  Real J = p[5] "Spin angular momentum about +z, in kg m^2/s times G/c^2 if geometric";
  Real r = sqrt(x[2]^2 + x[3]^2 + x[4]^2);
  Real phi = GeneralRelativity.Perturbation.Potentials.oblateJ2(x[2:4], {GM, R, J2});
  Real h0i[3];
algorithm
  // Lense-Thirring: h_0i = -2 G (J x r)_i / (c^3 r^3). With J given as G*J/c^2 in
  // geometric form this is -2 (J x r)_i / (c r^3); callers supply SI J and SI G.
  h0i := -2.0 * GeneralRelativity.Constants.G / (c^3 * r^3) * cross({0.0, 0.0, J}, x[2:4]);
  g := GeneralRelativity.Perturbation.weakFieldMetric(phi, h0i, c);
  annotation(Documentation(info="<html>
<p>Post-Newtonian Earth to first order: <code>g<sub>00</sub> = -(1 + 2&Phi;/c&sup2;)</code>,
<code>g<sub>ij</sub> = (1 - 2&Phi;/c&sup2;)&delta;<sub>ij</sub></code> with the zonal potential
<a href=\"modelica://GeneralRelativity.Perturbation.Potentials.oblateJ2\">oblateJ2</a>, and
the gravitomagnetic term <code>g<sub>0i</sub> = -2G(J &times; r)<sub>i</sub>/(c&sup3;r&sup3;)</code>
responsible for Lense-Thirring precession. Set <code>J = 0</code> to drop frame
dragging. The <code>z</code> axis is the spin axis. The gravitomagnetic term
uses the SI Newton constant from
<a href=\"modelica://GeneralRelativity.Constants\">Constants</a>; pass SI units.</p>
<p>Finite differences of this metric are below double-precision resolution at
Earth scales (<code>&Phi;/c&sup2; &asymp; 7&times;10<sup>-10</sup></code>). Integrate
orbits with
<a href=\"modelica://GeneralRelativity.Perturbation.postNewtonianAcceleration\">postNewtonianAcceleration</a>
and evaluate tides with
<a href=\"modelica://GeneralRelativity.Perturbation.weakFieldCurvature\">weakFieldCurvature</a>;
use this function to inspect metric components or feed a visualization.</p>
</html>"));
end earth;
