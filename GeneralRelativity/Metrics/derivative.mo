within GeneralRelativity.Metrics;
function derivative "Central-difference metric gradient dg[k, i, j] = d g_ij / d x^k"
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric function";
  input Real x[4] "Coordinates";
  input Real p[:] "Metric parameters";
  input Real step[4] "Finite-difference step per coordinate";
  output Real dg[4, 4, 4];
protected
  Real xPlus[4];
  Real xMinus[4];
  Real gPlus[4, 4];
  Real gMinus[4, 4];
algorithm
  for k in 1:4 loop
    xPlus := x;
    xMinus := x;
    xPlus[k] := x[k] + step[k];
    xMinus[k] := x[k] - step[k];
    // Separate calls: OpenModelica's C backend cannot combine two functional-input
    // calls inside one slice assignment.
    gPlus := metric(xPlus, p);
    gMinus := metric(xMinus, p);
    dg[k, :, :] := (gPlus - gMinus) / (2.0 * step[k]);
  end for;
end derivative;
