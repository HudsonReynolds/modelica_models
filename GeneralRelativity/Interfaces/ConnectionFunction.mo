within GeneralRelativity.Interfaces;
partial function ConnectionFunction "Affine connection coefficients Gamma^i_jk(x) of a chart, not necessarily from a metric"
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Connection parameter vector; layout is documented by each connection";
  output Real Gamma[4, 4, 4] "Gamma[i, j, k] = Gamma^i_jk; the antisymmetric part in j, k is the torsion";
  annotation(Documentation(info="<html>
<p>The covariant derivative of a vector is
<code>&nabla;<sub>j</sub>v<sup>i</sup> = &part;<sub>j</sub>v<sup>i</sup> + &Gamma;<sup>i</sup><sub>jk</sub>v<sup>k</sup></code>,
so the derivative index is the second array index, as in the Levi-Civita
case returned by <a href=\"modelica://GeneralRelativity.Tensors.christoffel\">Tensors.christoffel</a>.
Use this signature for connections that are not built from a Lorentzian
metric, such as the Newton-Cartan connection of
<a href=\"modelica://GeneralRelativity.NewtonCartan\">NewtonCartan</a> or a
metric-affine connection with torsion.</p>
</html>"));
end ConnectionFunction;
