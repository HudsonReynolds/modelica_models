within GeneralRelativity.Perturbation.Potentials;
function twoBody "Phi = -GM1/|x| - GM2/|x - x2|; p = {GM1, GM2, x2, y2, z2}"
  extends GeneralRelativity.Interfaces.PotentialFunction;
protected
  Real d[3] = x - p[3:5];
algorithm
  phi := -p[1] / sqrt(x * x) - p[2] / sqrt(d * d);
end twoBody;
