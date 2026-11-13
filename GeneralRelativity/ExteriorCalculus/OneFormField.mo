within GeneralRelativity.ExteriorCalculus;
partial function OneFormField "1-form field w_a(x) in coordinate components"
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Field parameter vector";
  output Real w[4] "Covariant components w_a";
end OneFormField;
