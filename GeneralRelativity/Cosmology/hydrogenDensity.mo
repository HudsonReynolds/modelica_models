within GeneralRelativity.Cosmology;
function hydrogenDensity "Number density of hydrogen nuclei today, 1/m^3"
  input Parameters par;
  output Real nH(unit = "1/m3");
protected
  constant Real mH = 1.6735575e-27 "Mass of a hydrogen atom";
  Real rhoCrit = 3.0 * (100.0e3 / GeneralRelativity.Constants.Mpc)^2 / (8.0 * GeneralRelativity.Constants.pi * GeneralRelativity.Constants.G) "Critical density / h^2";
algorithm
  nH := (1.0 - par.YHe) * par.omegaBh2 * rhoCrit / mH;
end hydrogenDensity;
