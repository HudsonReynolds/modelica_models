within GeneralRelativity.ExteriorCalculus;
partial function ThreeFormField "3-form field H_abc(x) as a totally antisymmetric component array"
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Field parameter vector";
  output Real H[4, 4, 4] "Totally antisymmetric components H_abc, H = (1/3!) H_abc dx^a ^ dx^b ^ dx^c";
end ThreeFormField;
