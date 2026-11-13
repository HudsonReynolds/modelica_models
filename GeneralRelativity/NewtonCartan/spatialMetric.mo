within GeneralRelativity.NewtonCartan;
function spatialMetric "Euclidean metric gamma_ij of an instant, the inverse of the spatial block of h^ab"
  input Real x[4] "Coordinates {ct, x, y, z} or {ct, r, theta, phi}";
  input Boolean spherical = false "true for {ct, r, theta, phi}";
  output Real gamma[3, 3] "gamma_ij for the spatial coordinates x[2:4]";
algorithm
  gamma := if spherical then diagonal({1.0, x[2]^2, (x[2] * sin(x[3]))^2}) else identity(3);
end spatialMetric;
