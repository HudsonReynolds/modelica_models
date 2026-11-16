within GeneralRelativity.Lensing;
function angularDiameterDistance "Angular-diameter distance from redshift z1 to z2 in a matter + Lambda Friedmann-Lemaitre universe"
  input Real z1 = 0.0 "Redshift of the nearer object (0 for the observer)";
  input Real z2 "Redshift of the farther object, z2 >= z1";
  input Real H0 = 70.0 "Hubble constant in km/s/Mpc";
  input Real OmegaM = 0.3 "Matter density parameter";
  input Real OmegaLambda = 1.0 - OmegaM "Dark-energy density parameter; the default makes the universe flat";
  input Integer n = 2000 "Simpson intervals for the comoving-distance integral (even)";
  output Real D(unit = "m") "Angular-diameter distance D_A(z1, z2)";
protected
  constant Real c = GeneralRelativity.Constants.c;
  Real DH = c / (H0 * 1000.0 / GeneralRelativity.Constants.Mpc) "Hubble distance c/H0";
  Real OmegaK = 1.0 - OmegaM - OmegaLambda;
  Real chi[2] "Line-of-sight comoving distances to z1 and z2";
  Real DM[2] "Transverse comoving distances";
  Integer m = 2 * max(1, div(n, 2));
  Real h;
  Real sum;
  Real zi;
  Real z;
algorithm
  assert(z2 >= z1, "angularDiameterDistance: z2 must not be smaller than z1");
  for j in 1:2 loop
    z := if j == 1 then z1 else z2;
    h := z / m;
    sum := 0.0;
    for i in 0:m loop
      zi := i * h;
      sum := sum + (if i == 0 or i == m then 1.0 else if mod(i, 2) == 1 then 4.0 else 2.0)
        / sqrt(OmegaM * (1.0 + zi)^3 + OmegaK * (1.0 + zi)^2 + OmegaLambda);
    end for;
    chi[j] := DH * sum * h / 3.0;
    DM[j] := if OmegaK > 1.0e-12 then DH / sqrt(OmegaK) * sinh(sqrt(OmegaK) * chi[j] / DH)
      else if OmegaK < -1.0e-12 then DH / sqrt(-OmegaK) * sin(sqrt(-OmegaK) * chi[j] / DH)
      else chi[j];
  end for;
  // Hogg (1999), eq. 19: exact for any curvature (reduces to (chi2 - chi1)/(1 + z2) when flat).
  D := (DM[2] * sqrt(1.0 + OmegaK * DM[1]^2 / DH^2) - DM[1] * sqrt(1.0 + OmegaK * DM[2]^2 / DH^2)) / (1.0 + z2);
  annotation(Documentation(info="<html>
<p>Line-of-sight comoving distance <code>&chi;(z) = D<sub>H</sub> &int;<sub>0</sub><sup>z</sup> dz'/E(z')</code> with
<code>E(z) = &radic;(&Omega;<sub>m</sub>(1+z)&sup3; + &Omega;<sub>k</sub>(1+z)&sup2; + &Omega;<sub>&Lambda;</sub>)</code>
(Simpson's rule), transverse comoving distance <code>D<sub>M</sub></code> for the curvature
<code>&Omega;<sub>k</sub></code>, and Hogg's formula for the distance between two redshifts. Radiation is
neglected. For <code>H<sub>0</sub> = 70</code>, <code>&Omega;<sub>m</sub> = 0.3</code>, flat:
<code>D<sub>A</sub>(0, 1) = 1651.9 Mpc</code>.</p>
</html>"));
end angularDiameterDistance;
