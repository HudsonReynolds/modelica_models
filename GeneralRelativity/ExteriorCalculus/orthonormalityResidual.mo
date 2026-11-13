within GeneralRelativity.ExteriorCalculus;
function orthonormalityResidual "Largest deviation of g(e_a, e_b) from eta_ab = diag(-1, 1, 1, 1)"
  input Real g[4, 4] "Covariant metric";
  input Real e[4, 4] "Tetrad e[a, mu] = e_a^mu";
  output Real residual "max |e g e^T - eta|";
protected
  Real gram[4, 4];
algorithm
  gram := e * g * transpose(e) - diagonal({-1.0, 1.0, 1.0, 1.0});
  residual := max(abs(gram[a, b]) for a in 1:4, b in 1:4);
end orthonormalityResidual;
