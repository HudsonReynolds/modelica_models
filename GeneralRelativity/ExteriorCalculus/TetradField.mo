within GeneralRelativity.ExteriorCalculus;
partial function TetradField "Orthonormal frame field: e[a, mu] = contravariant components of frame vector a"
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Tetrad parameter vector; usually the metric parameters";
  output Real e[4, 4] "e[a, mu] = e_a^mu with g(e_a, e_b) = eta_ab = diag(-1, 1, 1, 1); row 1 is timelike";
end TetradField;
