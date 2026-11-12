within GeneralRelativity.Perturbation.Potentials;
function pointMass "Phi = -GM / r; p = {GM}"
  extends GeneralRelativity.Interfaces.PotentialFunction;
algorithm
  phi := -p[1] / sqrt(x * x);
end pointMass;
