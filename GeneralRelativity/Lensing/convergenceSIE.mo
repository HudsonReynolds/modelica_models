within GeneralRelativity.Lensing;
function convergenceSIE "Convergence kappa = q b / (2 sqrt(q^2 x^2 + y^2)) of a singular isothermal ellipse"
  input Real theta[2](each unit = "rad") "Position in the lens frame";
  input Real b(unit = "rad") "Deflection scale";
  input Real q = 1.0 "Axis ratio";
  output Real kappa "Surface density over the critical density";
algorithm
  kappa := q * b / (2.0 * sqrt(q^2 * theta[1]^2 + theta[2]^2));
  annotation(Inline = true);
end convergenceSIE;
