within GeneralRelativity.Cosmology;
function heliumElectrons "Free electrons per hydrogen atom supplied by helium in Saha equilibrium"
  input Real z;
  input Real ne "Free-electron density (all species), 1/m^3";
  input Parameters par;
  output Real xHe;
protected
  constant Real kB = 1.380649e-23;
  constant Real hP = 6.62607015e-34;
  constant Real me = 9.1093837015e-31;
  constant Real eV = 1.602176634e-19;
  Real T = par.TCMB * (1.0 + z);
  Real fHe = par.YHe / (3.9715 * (1.0 - par.YHe)) "n_He / n_H (m_He = 3.9715 m_H)";
  Real th = (2.0 * GeneralRelativity.Constants.pi * me * kB * T / hP^2)^1.5 / max(ne, 1.0e-30);
  Real r1 = 4.0 * th * exp(-24.587387 * eV / (kB * T)) "n(HeII) / n(HeI)";
  Real r2 = th * exp(-54.417760 * eV / (kB * T)) "n(HeIII) / n(HeII)";
algorithm
  xHe := fHe * (r1 + 2.0 * r1 * r2) / (1.0 + r1 + r1 * r2);
end heliumElectrons;
