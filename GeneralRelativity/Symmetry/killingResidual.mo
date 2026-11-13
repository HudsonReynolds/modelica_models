within GeneralRelativity.Symmetry;
function killingResidual "Relative size of L_xi g: about the finite-difference error for a Killing vector, order one otherwise"
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric function";
  input GeneralRelativity.Interfaces.VectorFieldFunction field "Candidate Killing vector";
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Metric parameters";
  input Real q[:] "Field parameters";
  input Real step[4] "Finite-difference step per coordinate";
  output Real residual "max |(L_xi g)_ab| / max_ab sum of term magnitudes; 0 when every term vanishes";
protected
  Real L[4, 4];
  Real scale;
  Real largest;
algorithm
  (L, scale) := lieDerivativeMetric(metric, field, x, p, q, step);
  largest := max(abs(L[a, b]) for a in 1:4, b in 1:4);
  residual := if scale > 0.0 then largest / scale else 0.0;
  annotation(Documentation(info="<html>
<p>The ratio lies in <code>[0, 1]</code>. It is independent of the normalization
of <code>&xi;</code> and of the units of the chart, so one tolerance (for
example <code>1e-5</code> with step <code>1e-3</code>) serves all fields.
A field whose terms all vanish identically at <code>x</code>, such as
<code>&part;<sub>t</sub></code> on a static metric, gives exactly zero.</p>
</html>"));
end killingResidual;
