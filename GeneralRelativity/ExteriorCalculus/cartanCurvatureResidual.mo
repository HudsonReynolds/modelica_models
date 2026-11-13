within GeneralRelativity.ExteriorCalculus;
function cartanCurvatureResidual "Difference between Cartan's curvature 2-forms and the tetrad-projected coordinate Riemann tensor"
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric function";
  input TetradField tetrad "Orthonormal tetrad field for this metric";
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Metric parameters";
  input Real q[:] "Tetrad parameters";
  input Real step[4] "Finite-difference step per coordinate";
  output Real residual[4, 4, 4, 4] "curvatureForm - frameRiemann";
  output Real relative "max |residual| / max |frameRiemann|";
protected
  GeneralRelativity.Tensors.Curvature k;
  Real e[4, 4];
  Real projected[4, 4, 4, 4];
  Real Omega[4, 4, 4, 4];
  Real largest;
algorithm
  k := GeneralRelativity.Metrics.curvatureAt(metric, x, p, step);
  e := tetrad(x, q);
  projected := frameRiemann(k.riemann, e);
  Omega := curvatureForm(metric, tetrad, x, p, q, step);
  residual := Omega - projected;
  largest := max(abs(projected[a, b, c, d]) for a in 1:4, b in 1:4, c in 1:4, d in 1:4);
  relative := if largest > 0.0 then
    max(abs(residual[a, b, c, d]) for a in 1:4, b in 1:4, c in 1:4, d in 1:4) / largest
    else max(abs(residual[a, b, c, d]) for a in 1:4, b in 1:4, c in 1:4, d in 1:4);
  annotation(Documentation(info="<html>
<p>Checks <code>d&omega;<sup>a</sup><sub>b</sub> + &omega;<sup>a</sup><sub>c</sub>&and;&omega;<sup>c</sup><sub>b</sub>
= &frac12;R<sup>a</sup><sub>b&mu;&nu;</sub>dx<sup>&mu;</sup>&and;dx<sup>&nu;</sup></code>
against <a href=\"modelica://GeneralRelativity.Metrics.curvatureAt\">Metrics.curvatureAt</a>.
For flat space (<code>largest = 0</code>) the absolute residual is returned in
<code>relative</code>.</p>
</html>"));
end cartanCurvatureResidual;
