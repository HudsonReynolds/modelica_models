within GeneralRelativity.Cosmology;
function modeDerivatives "Right-hand side of the linear perturbation equations for one wavenumber (conformal Newtonian gauge)"
  input Real y[35] "a, phi, delta_c, theta_c, delta_b, theta_b, delta_g, theta_g, delta_nu, theta_nu, F_nu2..F_nu6, F_g2..F_g10, G_g0..G_g10";
  input Real eta(unit = "Mpc") "Conformal time";
  input Real k(unit = "1/Mpc") "Wavenumber";
  input Integer phase "0: photons and baryons tightly coupled; 1: photon Boltzmann hierarchy with Thomson scattering; 2: radiation perturbations dropped";
  input Real dtau(unit = "1/Mpc") "Thomson scattering rate n_e sigma_T a";
  input Real H0 "1/Mpc";
  input Real omegaGamma;
  input Real omegaNu;
  input Real omegaB;
  input Real omegaC;
  input Real omegaK;
  input Real omegaLambda;
  output Real dy[35];
protected
  constant Integer lg = 10 "Highest photon multipole";
  Real a = y[1];
  Real aDot = a2H(a, H0, omegaGamma + omegaNu, omegaB + omegaC, omegaK, omegaLambda);
  Real calH = aDot / a "Conformal Hubble rate a'/a";
  Boolean radiation = phase < 2;
  Real psi;
  Real drho "4 pi G a^2 delta rho";
  Real dphi;
  Real R = 0.75 * omegaB / omegaGamma * a "3 rho_b / (4 rho_gamma)";
  Real thetaDot;
  Real k2 = k * k;
  Real Pi "Photon quadrupole plus polarization, F_g2 + G_g0 + G_g2";
algorithm
  // Anisotropic stress of neutrinos (and of photons once they leave tight coupling):
  // k^2 (phi - psi) = 12 pi G a^2 sum (rho + P) sigma, sigma = F_2 / 2.
  psi := y[2] - (if radiation then 3.0 * H0^2 * (omegaNu * y[11] + (if phase == 1 then omegaGamma * y[16] else 0.0)) / (a^2 * k2) else 0.0);
  drho := 1.5 * H0^2 * ((omegaC * y[3] + omegaB * y[5]) / a + (if radiation then (omegaGamma * y[7] + omegaNu * y[9]) / a^2 else 0.0));
  if k * eta < 10.0 then
    // Outside and at the horizon: the 00 Einstein equation (Ma & Bertschinger eq. 23a),
    // k^2 phi + 3 H (phi' + H psi) = -4 pi G a^2 delta rho. (The 0i equation alone admits phi = const on
    // superhorizon scales and so misses the decay of phi across the radiation-matter transition.)
    dphi := -calH * psi - (k2 * y[2] + drho) / (3.0 * calH);
  else
    // Well inside the horizon the 00 equation pulls phi to the Poisson constraint at the rate k^2 / (3 H), which
    // is stiff; the 0i equation (MB eq. 23b), k^2 (phi' + H psi) = 4 pi G a^2 (rho + P) theta, is not.
    dphi := 1.5 * H0^2 * ((omegaC * y[4] + omegaB * y[6]) / a + (if radiation then 4.0 / 3.0 * (omegaGamma * y[8] + omegaNu * y[10]) / a^2 else 0.0)) / k2
      - calH * psi;
  end if;
  dy := zeros(35);
  dy[1] := aDot;
  dy[2] := dphi;
  dy[3] := -y[4] + 3.0 * dphi;
  dy[4] := -calH * y[4] + k2 * psi;
  dy[5] := -y[6] + 3.0 * dphi;
  if phase == 0 then
    // Tight coupling: theta_b = theta_gamma, a single fluid loaded by the baryons; no photon anisotropy.
    thetaDot := (-calH * R * y[6] + k2 * y[7] / 4.0) / (1.0 + R) + k2 * psi;
    dy[6] := thetaDot;
    dy[8] := thetaDot;
  elseif phase == 1 then
    // Baryons dragged by Thomson scattering (MB eq. 66), photons as a Boltzmann hierarchy (MB eq. 63-64).
    dy[6] := -calH * y[6] + k2 * psi + dtau / R * (y[8] - y[6]);
    dy[8] := k2 * (y[7] / 4.0 - y[16] / 2.0) + k2 * psi + dtau * (y[6] - y[8]);
    Pi := y[16] + y[25] + y[27];
    dy[16] := 8.0 / 15.0 * y[8] - 3.0 / 5.0 * k * y[17] - 0.9 * dtau * y[16] + 0.1 * dtau * (y[25] + y[27]);
    for l in 3:lg - 1 loop
      dy[l + 14] := k / (2 * l + 1) * (l * y[l + 13] - (l + 1) * y[l + 15]) - dtau * y[l + 14];
    end for;
    dy[lg + 14] := k * y[lg + 13] - (lg + 1) / eta * y[lg + 14] - dtau * y[lg + 14];
    // Polarization G_l at y[l + 25].
    dy[25] := -k * y[26] + dtau * (-y[25] + Pi / 2.0);
    dy[26] := k / 3.0 * (y[25] - 2.0 * y[27]) - dtau * y[26];
    dy[27] := k / 5.0 * (2.0 * y[26] - 3.0 * y[28]) + dtau * (-y[27] + Pi / 10.0);
    for l in 3:lg - 1 loop
      dy[l + 25] := k / (2 * l + 1) * (l * y[l + 24] - (l + 1) * y[l + 26]) - dtau * y[l + 25];
    end for;
    dy[lg + 25] := k * y[lg + 24] - (lg + 1) / eta * y[lg + 25] - dtau * y[lg + 25];
  else
    dy[6] := -calH * y[6] + k2 * psi;
  end if;
  if radiation then
    dy[7] := -4.0 / 3.0 * y[8] + 4.0 * dphi;
    dy[9] := -4.0 / 3.0 * y[10] + 4.0 * dphi;
    dy[10] := k2 * (y[9] / 4.0 - y[11] / 2.0) + k2 * psi;
    dy[11] := 8.0 / 15.0 * y[10] - 3.0 / 5.0 * k * y[12];
    for l in 3:5 loop
      dy[l + 9] := k / (2 * l + 1) * (l * y[l + 8] - (l + 1) * y[l + 10]);
    end for;
    dy[15] := k * y[14] - 7.0 / eta * y[15];
  end if;
  annotation(Documentation(info="<html>
<p>Ma &amp; Bertschinger (ApJ 455, 7, 1995), conformal Newtonian gauge
<code>ds&sup2; = a&sup2;[&minus;(1 + 2&psi;)d&eta;&sup2; + (1 &minus; 2&phi;)dx&sup2;]</code>:</p>
<ul>
<li>cold dark matter: <code>&delta;<sub>c</sub>' = &minus;&theta;<sub>c</sub> + 3&phi;'</code>, <code>&theta;<sub>c</sub>' = &minus;(a'/a)&theta;<sub>c</sub> + k&sup2;&psi;</code>;</li>
<li>photons and baryons, while Thomson scattering is fast (phase 0), as one fluid:
<code>&theta;' = [&minus;(a'/a)R&theta; + k&sup2;&delta;<sub>&gamma;</sub>/4]/(1 + R) + k&sup2;&psi;</code> &mdash; an oscillator with
sound speed <code>1/&radic;(3(1 + R))</code> driven by the potentials: the acoustic waves;</li>
<li>then (phase 1) photons as a Boltzmann hierarchy to <code>l = 10</code> with polarization, scattering off the
baryons at the rate <code>&tau;'</code>: diffusion (Silk) damping and decoupling follow from the equations;</li>
<li>massless neutrinos: hierarchy to <code>l = 6</code> with the free-streaming truncation; neutrino (and photon)
anisotropic stress separates <code>&phi;</code> and <code>&psi;</code>;</li>
<li><code>&phi;</code> from the 00 Einstein equation outside the horizon and the 0i equation inside it (where the
00 equation is stiff).</li>
</ul>
<p>Once a mode is deep inside the horizon after recombination (phase 2) the radiation perturbations are
dropped, the radiation-streaming approximation of CAMB and CLASS, so the integration to today is cheap.</p>
</html>"));
end modeDerivatives;
