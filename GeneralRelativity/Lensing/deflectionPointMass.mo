within GeneralRelativity.Lensing;
function deflectionPointMass "Scaled deflection thetaE^2 theta / |theta|^2 of a point mass at the origin"
  input Real theta[2](each unit = "rad") "Image-plane position relative to the lens";
  input Real thetaE(unit = "rad") "Einstein radius";
  output Real alpha[2](each unit = "rad");
protected
  Real r2 = theta * theta;
algorithm
  alpha := if r2 > 0.0 then thetaE^2 * theta / r2 else zeros(2);
end deflectionPointMass;
