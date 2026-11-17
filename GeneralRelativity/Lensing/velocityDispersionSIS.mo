within GeneralRelativity.Lensing;
function velocityDispersionSIS "Velocity dispersion c sqrt(thetaE Ds / (4 pi Dls)) of the isothermal sphere with a given Einstein radius: the inverse of einsteinRadiusSIS"
  input Real thetaE(unit = "rad") "Einstein radius";
  input Real Ds(unit = "m") "Distance to the source";
  input Real Dls(unit = "m") "Distance from the lens to the source";
  output Real sigma(unit = "m/s") "Line-of-sight velocity dispersion";
algorithm
  sigma := GeneralRelativity.Constants.c * sqrt(thetaE * Ds / (4.0 * GeneralRelativity.Constants.pi * Dls));
  annotation(Inline = true);
end velocityDispersionSIS;
