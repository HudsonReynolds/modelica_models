within GeneralRelativity.Metrics;
function ppn "Parametrized post-Newtonian point-mass metric in isotropic Cartesian coordinates {ct, x, y, z}; p = {rs, gamma, beta}"
  extends GeneralRelativity.Interfaces.MetricFunction;
protected
  Real rs = p[1] "Schwarzschild radius 2GM/c^2";
  Real gamma = p[2] "PPN gamma: space curvature per unit mass (1 in GR)";
  Real beta = p[3] "PPN beta: nonlinearity of superposition (1 in GR)";
  Real U = rs / (2.0 * sqrt(x[2]^2 + x[3]^2 + x[4]^2)) "Dimensionless Newtonian potential GM/(c^2 r)";
  Real gs = 1.0 + 2.0 * gamma * U;
algorithm
  g := diagonal({-1.0 + 2.0 * U - 2.0 * beta * U^2, gs, gs, gs});
  annotation(Documentation(info="<html>
<p>The static, spherically symmetric part of the PPN metric in the standard
isotropic gauge (Will, <i>Theory and Experiment in Gravitational Physics</i>):
<code>g<sub>00</sub> = -1 + 2U - 2&beta;U&sup2;</code>,
<code>g<sub>ij</sub> = (1 + 2&gamma;U)&delta;<sub>ij</sub></code>, <code>g<sub>0i</sub> = 0</code>,
with <code>U = GM/(c&sup2;r) = r<sub>s</sub>/(2r)</code>. The metric is
accurate to <code>O(U&sup2;)</code> in <code>g<sub>00</sub></code> and <code>O(U)</code> in
<code>g<sub>ij</sub></code>; with <code>&gamma; = &beta; = 1</code> it agrees with
<a href=\"modelica://GeneralRelativity.Metrics.isotropicSchwarzschild\">isotropicSchwarzschild</a>
to that order (the remainders are <code>3U&sup3;/2</code> and <code>3U&sup2;/2</code>). Brans-Dicke
theory has <code>&gamma; = (1 + &omega;)/(2 + &omega;)</code>, <code>&beta; = 1</code>.
The classical tests are in
<a href=\"modelica://GeneralRelativity.Perturbation\">Perturbation</a>
(<code>ppnPerihelionAdvance</code>, <code>ppnLightDeflection</code>, <code>ppnShapiroDelay</code>).</p>
</html>"));
end ppn;
