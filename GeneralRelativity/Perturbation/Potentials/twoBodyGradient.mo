within GeneralRelativity.Perturbation.Potentials;
function twoBodyGradient "Gradient of twoBody; p = {GM1, GM2, x2, y2, z2}"
  extends GeneralRelativity.Interfaces.PotentialGradientFunction;
protected
  Real d[3] = x - p[3:5];
  Real r1 = sqrt(x * x);
  Real r2 = sqrt(d * d);
algorithm
  gradient := p[1] * x / r1^3 + p[2] * d / r2^3;
end twoBodyGradient;
