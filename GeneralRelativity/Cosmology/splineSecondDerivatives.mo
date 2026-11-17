within GeneralRelativity.Cosmology;
function splineSecondDerivatives "Second derivatives of natural cubic splines through each column of f(x)"
  input Real x[:] "Increasing nodes";
  input Real f[size(x, 1), :] "Values, one spline per column";
  output Real M[size(x, 1), size(f, 2)];
protected
  Integer n = size(x, 1);
  Real c[size(x, 1)];
  Real d[size(x, 1), size(f, 2)];
  Real w;
algorithm
  // Thomas algorithm for (h_{i-1} M_{i-1} + 2 (h_{i-1} + h_i) M_i + h_i M_{i+1}) / 6 = slope differences,
  // with M_1 = M_n = 0.
  M := zeros(n, size(f, 2));
  c := zeros(n);
  d := zeros(n, size(f, 2));
  for i in 2:n - 1 loop
    w := 2.0 * (x[i + 1] - x[i - 1]) - (x[i] - x[i - 1]) * c[i - 1];
    c[i] := (x[i + 1] - x[i]) / w;
    for col in 1:size(f, 2) loop
      d[i, col] := (6.0 * ((f[i + 1, col] - f[i, col]) / (x[i + 1] - x[i]) - (f[i, col] - f[i - 1, col]) / (x[i] - x[i - 1]))
        - (x[i] - x[i - 1]) * d[i - 1, col]) / w;
    end for;
  end for;
  for i in n - 1:-1:2 loop
    for col in 1:size(f, 2) loop
      M[i, col] := d[i, col] - c[i] * M[i + 1, col];
    end for;
  end for;
end splineSecondDerivatives;
