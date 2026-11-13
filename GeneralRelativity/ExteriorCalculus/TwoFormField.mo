within GeneralRelativity.ExteriorCalculus;
partial function TwoFormField "2-form field F_ab(x) as an antisymmetric component array"
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Field parameter vector";
  output Real F[4, 4] "Antisymmetric components F_ab, F = (1/2) F_ab dx^a ^ dx^b";
end TwoFormField;
