within GeneralRelativity.Metrics;
function secondDerivative "Central-difference metric Hessian ddg[k, l, i, j] = d^2 g_ij / (d x^k d x^l)"
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric function";
  input Real x[4] "Coordinates";
  input Real p[:] "Metric parameters";
  input Real step[4] "Finite-difference step per coordinate";
  output Real ddg[4, 4, 4, 4];
protected
  Real g0[4, 4] = metric(x, p);
  Real xa[4];
  Real xb[4];
  Real xc[4];
  Real xd[4];
  Real ga[4, 4];
  Real gb[4, 4];
  Real gc[4, 4];
  Real gd[4, 4];
  Real mixed[4, 4];
algorithm
  for k in 1:4 loop
    xa := x;
    xb := x;
    xa[k] := x[k] + step[k];
    xb[k] := x[k] - step[k];
    ga := metric(xa, p);
    gb := metric(xb, p);
    ddg[k, k, :, :] := (ga - 2.0 * g0 + gb) / step[k]^2;
  end for;
  for k in 1:4 loop
    for l in k + 1:4 loop
      xa := x;
      xb := x;
      xc := x;
      xd := x;
      xa[k] := x[k] + step[k];
      xa[l] := x[l] + step[l];
      xb[k] := x[k] + step[k];
      xb[l] := x[l] - step[l];
      xc[k] := x[k] - step[k];
      xc[l] := x[l] + step[l];
      xd[k] := x[k] - step[k];
      xd[l] := x[l] - step[l];
      ga := metric(xa, p);
      gb := metric(xb, p);
      gc := metric(xc, p);
      gd := metric(xd, p);
      mixed := (ga - gb - gc + gd) / (4.0 * step[k] * step[l]);
      ddg[k, l, :, :] := mixed;
      ddg[l, k, :, :] := mixed;
    end for;
  end for;
end secondDerivative;
