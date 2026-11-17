within GeneralRelativity.Cosmology;
function splineEvaluate "Natural cubic spline through (x, f) with second derivatives M, evaluated at xq"
  input Real x[:] "Increasing nodes";
  input Real f[size(x, 1)];
  input Real M[size(x, 1)] "From splineSecondDerivatives";
  input Real xq;
  output Real y;
protected
  Integer n = size(x, 1);
  Integer lo = 1;
  Integer hi = size(x, 1);
  Integer mid;
  Real h;
  Real A;
  Real B;
algorithm
  if xq <= x[1] then
    y := f[1];
  elseif xq >= x[n] then
    y := f[n];
  else
    while hi - lo > 1 loop
      mid := div(lo + hi, 2);
      if x[mid] <= xq then
        lo := mid;
      else
        hi := mid;
      end if;
    end while;
    h := x[hi] - x[lo];
    A := (x[hi] - xq) / h;
    B := 1.0 - A;
    y := A * f[lo] + B * f[hi] + ((A^3 - A) * M[lo] + (B^3 - B) * M[hi]) * h^2 / 6.0;
  end if;
end splineEvaluate;
