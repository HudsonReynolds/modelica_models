within GeneralRelativity.Interfaces;
partial function MetricFunction "Covariant metric g_ij(x) in the (-+++) signature"
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Metric parameter vector; layout is documented by each metric";
  output Real g[4, 4] "Covariant metric components g_ij";
end MetricFunction;
