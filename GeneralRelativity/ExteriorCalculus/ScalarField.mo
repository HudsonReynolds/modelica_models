within GeneralRelativity.ExteriorCalculus;
partial function ScalarField "Scalar field (0-form) f(x)"
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Field parameter vector";
  output Real f "Value of the 0-form";
end ScalarField;
