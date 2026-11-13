within GeneralRelativity.Symmetry;
function conformalKillingResidual "Relative size of L_xi g - (2/n)(div xi) g with n = 4; small for a conformal Killing vector"
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric function";
  input GeneralRelativity.Interfaces.VectorFieldFunction field "Candidate conformal Killing vector";
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Metric parameters";
  input Real q[:] "Field parameters";
  input Real step[4] "Finite-difference step per coordinate";
  output Real residual "max |L_xi g - (div xi / 2) g| / (Lie-derivative term scale + max |(div xi / 2) g|)";
protected
  Real g[4, 4];
  Real L[4, 4];
  Real scale;
  Real div;
  Real trace[4, 4];
algorithm
  g := metric(x, p);
  (L, scale) := lieDerivativeMetric(metric, field, x, p, q, step);
  div := divergence(metric, field, x, p, q, step);
  trace := (2.0 / 4.0) * div * g;
  scale := scale + max(abs(trace[a, b]) for a in 1:4, b in 1:4);
  residual := if scale > 0.0 then max(abs(L[a, b] - trace[a, b]) for a in 1:4, b in 1:4) / scale else 0.0;
  annotation(Documentation(info="<html>
<p>The trace-free part of <code>&Lscr;<sub>&xi;</sub>g</code>. Killing vectors
are conformal Killing vectors with zero divergence; the
<a href=\"modelica://GeneralRelativity.Symmetry.dilation\">dilation</a>
<code>x<sup>a</sup>&part;<sub>a</sub></code> of Minkowski space is conformal Killing
(<code>&Lscr;<sub>&xi;</sub>&eta; = 2&eta;</code>) but not Killing.</p>
</html>"));
end conformalKillingResidual;
