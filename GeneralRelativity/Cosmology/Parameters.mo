within GeneralRelativity.Cosmology;
record Parameters "Flat or curved Lambda-CDM parameters; defaults are the Planck 2018 best fit (TT,TE,EE+lowE+lensing)"
  Real h = 0.6736 "Hubble constant H0 / (100 km/s/Mpc)";
  Real omegaBh2 = 0.02237 "Baryon density Omega_b h^2";
  Real omegaCh2 = 0.1200 "Cold dark matter density Omega_c h^2";
  Real omegaK = 0.0 "Spatial curvature Omega_k (Lambda follows from closure)";
  Real TCMB(unit = "K") = 2.7255 "CMB temperature today";
  Real Neff = 3.046 "Effective number of (massless) neutrino species";
  Real YHe = 0.2454 "Primordial helium mass fraction";
  Real nS = 0.9649 "Scalar spectral index";
  Real AS = 2.100e-9 "Amplitude of the primordial curvature power spectrum at the pivot";
  Real kPivot(unit = "1/Mpc") = 0.05 "Pivot scale of AS and nS";
  Real tauReio = 0.0544 "Optical depth to reionization (sets the reionization redshift)";
  annotation(Documentation(info="<html>
<p>Planck Collaboration, <i>Planck 2018 results. VI. Cosmological parameters</i>, A&amp;A 641, A6 (2020),
table 2 (TT,TE,EE+lowE+lensing). Neutrinos are treated as massless here (Planck's base model has one
massive species with 0.06 eV), which lengthens the age by about 0.02 Gyr.</p>
</html>"));
end Parameters;
