within GeneralRelativity.Cosmology;
function reionizationFraction "Free electrons per hydrogen atom from reionization (CAMB's tanh model, helium doubly ionized at z = 3.5)"
  input Real z;
  input Real zReio "Mid-point of hydrogen reionization";
  input Parameters par;
  output Real x;
protected
  Real fHe = par.YHe / (3.9715 * (1.0 - par.YHe));
  Real y = (1.0 + z)^1.5;
  Real yRe = (1.0 + zReio)^1.5;
  Real dy = 1.5 * sqrt(1.0 + zReio) * 0.5 "Width dz = 0.5";
algorithm
  x := (1.0 + fHe) / 2.0 * (1.0 + tanh((yRe - y) / dy)) + fHe / 2.0 * (1.0 + tanh((3.5 - z) / 0.5));
end reionizationFraction;
