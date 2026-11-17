within GeneralRelativity.Cosmology;
function lastScatteringRedshift "Redshift of the peak of the visibility function"
  input Parameters par;
  output Real zStar;
protected
  Real z[6000], a[6000], eta[6000], xe[6000], dtau[6000], tau[6000], g[6000], kD[6000];
  Real etaStar;
  Real zReio;
algorithm
  (z, a, eta, xe, dtau, tau, g, kD, zStar, etaStar, zReio) := thermalHistory(par, 6000);
end lastScatteringRedshift;
