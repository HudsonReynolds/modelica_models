within GeneralRelativity.Perturbation.Potentials;
function pointMassGradient "grad Phi = GM x / r^3; p = {GM}"
  extends GeneralRelativity.Interfaces.PotentialGradientFunction;
protected
  Real r = sqrt(x * x);
algorithm
  gradient := p[1] * x / r^3;
end pointMassGradient;
