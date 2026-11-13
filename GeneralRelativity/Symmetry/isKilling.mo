within GeneralRelativity.Symmetry;
function isKilling "True when killingResidual is below a tolerance at the given event"
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric function";
  input GeneralRelativity.Interfaces.VectorFieldFunction field "Candidate Killing vector";
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Metric parameters";
  input Real q[:] "Field parameters";
  input Real step[4] "Finite-difference step per coordinate";
  input Real tolerance = 1.0e-5 "Largest accepted relative residual";
  output Boolean result "killingResidual < tolerance";
algorithm
  result := killingResidual(metric, field, x, p, q, step) < tolerance;
  annotation(Documentation(info="<html>
<p>A pointwise test: a field can satisfy Killing's equation at isolated events by
accident, so check several events (as the tests do) before concluding that a
field is a symmetry.</p>
</html>"));
end isKilling;
