within GeneralRelativity.Cosmology;
function reionizationOpticalDepth "Thomson optical depth of the reionized electrons, int n_e sigma_T c dt"
  input Real zReio;
  input Parameters par;
  input Integer n = 4000 "Simpson intervals in z";
  output Real tau;
protected
  constant Real sigmaT = 6.6524587321e-29;
  Integer m = 2 * max(1, div(n, 2));
  Real zMax = zReio + 15.0;
  Real dz = zMax / m;
  Real z;
  Real s = 0.0;
algorithm
  for i in 0:m loop
    z := i * dz;
    s := s + (if i == 0 or i == m then 1.0 else if mod(i, 2) == 1 then 4.0 else 2.0)
      * reionizationFraction(z, zReio, par) * hydrogenDensity(par) * (1.0 + z)^2 * sigmaT * GeneralRelativity.Constants.c
      / (hubble(1.0 / (1.0 + z), par) * 1000.0 / GeneralRelativity.Constants.Mpc);
  end for;
  tau := s * dz / 3.0;
end reionizationOpticalDepth;
