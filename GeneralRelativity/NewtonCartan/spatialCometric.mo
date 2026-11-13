within GeneralRelativity.NewtonCartan;
function spatialCometric "Degenerate spatial cometric h^ab of Euclidean space in a Cartesian or spherical chart"
  input Real x[4] "Coordinates {ct, x, y, z} or {ct, r, theta, phi}";
  input Boolean spherical = false "true for {ct, r, theta, phi}";
  output Real h[4, 4] "h^ab; h^ab tau_b = 0";
algorithm
  h := zeros(4, 4);
  if spherical then
    h[2, 2] := 1.0;
    h[3, 3] := 1.0 / x[2]^2;
    h[4, 4] := 1.0 / (x[2] * sin(x[3]))^2;
  else
    h[2, 2] := 1.0;
    h[3, 3] := 1.0;
    h[4, 4] := 1.0;
  end if;
end spatialCometric;
