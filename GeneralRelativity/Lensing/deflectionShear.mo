within GeneralRelativity.Lensing;
function deflectionShear "Scaled deflection of an external shear gamma at position angle phi (the tidal field of the lens's environment)"
  input Real theta[2](each unit = "rad") "Position relative to the lens";
  input Real gamma "Shear strength";
  input Real phi(unit = "rad") = 0.0 "Shear position angle";
  output Real alpha[2](each unit = "rad");
protected
  Real g1 = gamma * cos(2.0 * phi);
  Real g2 = gamma * sin(2.0 * phi);
algorithm
  // Lensing potential psi = (g1 (x^2 - y^2) + 2 g2 x y) / 2.
  alpha := {g1 * theta[1] + g2 * theta[2], g2 * theta[1] - g1 * theta[2]};
  annotation(Inline = true);
end deflectionShear;
