within GeneralRelativity.ExteriorCalculus;
function exteriorDerivative0 "Exterior derivative of a scalar field, (df)_a = d_a f, by central differences"
  input ScalarField field "0-form field";
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Field parameters";
  input Real step[4] "Finite-difference step per coordinate";
  output Real w[4] "1-form df";
protected
  Real xPlus[4];
  Real xMinus[4];
  Real fPlus;
  Real fMinus;
algorithm
  for a in 1:4 loop
    xPlus := x;
    xMinus := x;
    xPlus[a] := x[a] + step[a];
    xMinus[a] := x[a] - step[a];
    fPlus := field(xPlus, p);
    fMinus := field(xMinus, p);
    w[a] := (fPlus - fMinus) / (2.0 * step[a]);
  end for;
end exteriorDerivative0;
