within GeneralRelativity.Perturbation.Potentials;
function oblateJ2 "Zonal potential Phi = -GM/r [1 - J2 (R/r)^2 P2(z/r)], spin axis +z; p = {GM, R, J2}"
  extends GeneralRelativity.Interfaces.PotentialFunction;
protected
  Real GM = p[1];
  Real R = p[2];
  Real J2 = p[3];
  Real r2 = x * x;
  Real r = sqrt(r2);
  Real sinLat2 = x[3]^2 / r2;
algorithm
  phi := -GM / r * (1.0 - J2 * (R / r)^2 * 0.5 * (3.0 * sinLat2 - 1.0));
end oblateJ2;
