within GeneralRelativity.Lensing;
function criticalDensity "Critical surface density c^2 Ds / (4 pi G Dl Dls); the convergence is Sigma / Sigma_cr"
  input Real Dl(unit = "m") "Angular-diameter distance to the lens";
  input Real Ds(unit = "m") "Angular-diameter distance to the source";
  input Real Dls(unit = "m") "Angular-diameter distance from the lens to the source";
  output Real SigmaCr(unit = "kg/m2");
algorithm
  SigmaCr := GeneralRelativity.Constants.c^2 * Ds / (4.0 * GeneralRelativity.Constants.pi * GeneralRelativity.Constants.G * Dl * Dls);
  annotation(Inline = true);
end criticalDensity;
