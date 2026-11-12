within GeneralRelativity.Metrics;
function minkowski "Flat metric diag(-1, 1, 1, 1) in Cartesian coordinates {ct, x, y, z}"
  extends GeneralRelativity.Interfaces.MetricFunction;
algorithm
  g := [-1.0, 0.0, 0.0, 0.0;
         0.0, 1.0, 0.0, 0.0;
         0.0, 0.0, 1.0, 0.0;
         0.0, 0.0, 0.0, 1.0];
end minkowski;
