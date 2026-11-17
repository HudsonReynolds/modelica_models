within GeneralRelativity.Cosmology;
function hydrogenSaha "Ionized fraction of hydrogen in Saha equilibrium (electrons from hydrogen only)"
  input Real z;
  input Parameters par;
  output Real x;
protected
  constant Real kB = 1.380649e-23;
  constant Real hP = 6.62607015e-34;
  constant Real me = 9.1093837015e-31;
  constant Real B1 = 13.605693 * 1.602176634e-19;
  Real T = par.TCMB * (1.0 + z);
  Real S = (2.0 * GeneralRelativity.Constants.pi * me * kB * T / hP^2)^1.5 * exp(-B1 / (kB * T)) / (hydrogenDensity(par) * (1.0 + z)^3);
algorithm
  // x^2 / (1 - x) = S
  x := 2.0 / (1.0 + sqrt(1.0 + 4.0 / S));
end hydrogenSaha;
