within GeneralRelativity.Metrics;
function schwarzschildChristoffel "Closed-form Christoffel symbols of the Schwarzschild metric, for reference and tests"
  input Real x[4] "Coordinates {ct, r, theta, phi}";
  input Real rs "Schwarzschild radius";
  output Real Gamma[4, 4, 4] "Gamma[i, j, k] = Gamma^i_jk";
protected
  Real r = x[2];
  Real theta = x[3];
  Real f = 1.0 - rs / r;
algorithm
  Gamma := zeros(4, 4, 4);
  Gamma[1, 1, 2] := rs / (2.0 * r^2 * f);
  Gamma[1, 2, 1] := Gamma[1, 1, 2];
  Gamma[2, 1, 1] := rs * f / (2.0 * r^2);
  Gamma[2, 2, 2] := -rs / (2.0 * r^2 * f);
  Gamma[2, 3, 3] := -r * f;
  Gamma[2, 4, 4] := -r * f * sin(theta)^2;
  Gamma[3, 2, 3] := 1.0 / r;
  Gamma[3, 3, 2] := 1.0 / r;
  Gamma[3, 4, 4] := -sin(theta) * cos(theta);
  Gamma[4, 2, 4] := 1.0 / r;
  Gamma[4, 4, 2] := 1.0 / r;
  Gamma[4, 3, 4] := cos(theta) / sin(theta);
  Gamma[4, 4, 3] := Gamma[4, 3, 4];
end schwarzschildChristoffel;
