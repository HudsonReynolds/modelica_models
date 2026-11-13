within GeneralRelativity.Metrics;
function isotropicSchwarzschild "Schwarzschild metric in isotropic Cartesian coordinates {ct, x, y, z}; p = {rs}"
  extends GeneralRelativity.Interfaces.MetricFunction;
protected
  Real rs = p[1] "Schwarzschild radius 2GM/c^2";
  Real rho = sqrt(x[2]^2 + x[3]^2 + x[4]^2) "Isotropic radius";
  Real q = rs / (4.0 * rho);
  Real psi4 = (1.0 + q)^4 "Conformal factor of the spatial metric";
algorithm
  g := diagonal({-((1.0 - q) / (1.0 + q))^2, psi4, psi4, psi4});
  annotation(Documentation(info="<html>
<p><code>ds&sup2; = -((1 - r<sub>s</sub>/4&rho;)/(1 + r<sub>s</sub>/4&rho;))&sup2; c&sup2;dt&sup2;
+ (1 + r<sub>s</sub>/4&rho;)<sup>4</sup>(dx&sup2; + dy&sup2; + dz&sup2;)</code> with
<code>&rho; = sqrt(x&sup2; + y&sup2; + z&sup2;)</code>. The areal (Schwarzschild) radius is
<code>r = &rho;(1 + r<sub>s</sub>/4&rho;)&sup2;</code>; the horizon is the sphere
<code>&rho; = r<sub>s</sub>/4</code>, and the chart covers the exterior twice
(<code>&rho; &gt; r<sub>s</sub>/4</code> and its inversion image).
The spatial metric is conformally flat, which makes this chart the strong-field
completion of the
<a href=\"modelica://GeneralRelativity.Metrics.ppn\">PPN</a> metric: expanding in
<code>U = r<sub>s</sub>/(2&rho;)</code> gives
<code>g<sub>00</sub> = -1 + 2U - 2U&sup2; + O(U&sup3;)</code> and
<code>g<sub>ij</sub> = (1 + 2U + O(U&sup2;))&delta;<sub>ij</sub></code>, i.e.
<code>&gamma; = &beta; = 1</code>.</p>
</html>"));
end isotropicSchwarzschild;
