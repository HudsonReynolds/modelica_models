within GeneralRelativity.Cosmology;
function thermalHistory "Ionization history, conformal time, optical depth, visibility function and diffusion-damping scale on a redshift grid"
  input Parameters par;
  input Integer n = 6000 "Grid points, uniform in ln(1 + z) from zMax down to 0";
  input Real zMax = 8000.0;
  output Real z[n] "Redshift (decreasing)";
  output Real a[n] "Scale factor";
  output Real eta[n](each unit = "Mpc") "Conformal time";
  output Real xe[n] "Free electrons per hydrogen atom (hydrogen, helium, reionization)";
  output Real dtau[n](each unit = "1/Mpc") "Thomson scattering rate per unit conformal time, n_e sigma_T a";
  output Real tau[n] "Optical depth from here to today";
  output Real g[n](each unit = "1/Mpc") "Visibility function dtau exp(-tau): where today's CMB photons last scattered";
  output Real kD[n](each unit = "1/Mpc") "Diffusion (Silk) damping wavenumber";
  output Real zStar "Redshift of the peak of the visibility function (last scattering)";
  output Real etaStar(unit = "Mpc") "Conformal time at zStar";
  output Real zReio "Reionization redshift";
protected
  constant Real sigmaT = 6.6524587321e-29;
  constant Real zSaha = 1600.0 "Hydrogen follows Saha above this redshift";
  Real H0, Og, On, Ob, Oc, OL;
  Real nH0;
  Real u;
  Real xH[n];
  Real xHe;
  Real zc;
  Real zn;
  Real dz;
  Real xc;
  Real fc;
  Real xn;
  Real fn;
  Real F;
  Real dF;
  Real eps;
  Integer m;
  Real R;
  Real Ri;
  Real Rj;
  Real wi;
  Real wj;
  Real sumD;
  Integer iMax;
  Real gl;
  Real gm;
  Real gr;
  Real shift;
algorithm
  (H0, Og, On, Ob, Oc, OL) := densities(par);
  nH0 := hydrogenDensity(par);
  zReio := reionizationRedshift(par);
  for i in 1:n loop
    u := log(1.0 + zMax) * (n - i) / (n - 1);
    z[i] := exp(u) - 1.0;
    a[i] := 1.0 / (1.0 + z[i]);
  end for;
  // Hydrogen: Saha while recombination is in equilibrium, then the Peebles equation (implicit trapezoid,
  // A-stable for the stiff early phase).
  for i in 1:n loop
    if z[i] >= zSaha then
      xH[i] := hydrogenSaha(z[i], par);
    else
      m := integer(ceil((z[i - 1] - z[i]) / 0.5));
      dz := (z[i] - z[i - 1]) / m;
      xc := xH[i - 1];
      zc := z[i - 1];
      for s in 1:m loop
        zn := zc + dz;
        fc := peeblesRate(zc, xc, heliumElectrons(zc, xc * nH0 * (1.0 + zc)^3, par), par);
        xn := xc + dz * fc;
        for it in 1:20 loop
          fn := peeblesRate(zn, xn, heliumElectrons(zn, xn * nH0 * (1.0 + zn)^3, par), par);
          F := xn - xc - 0.5 * dz * (fc + fn);
          eps := 1.0e-7 * max(xn, 1.0e-6);
          dF := 1.0 - 0.5 * dz * (peeblesRate(zn, xn + eps, heliumElectrons(zn, (xn + eps) * nH0 * (1.0 + zn)^3, par), par) - fn) / eps;
          xn := min(max(xn - F / dF, 1.0e-12), 1.0);
          if abs(F) < 1.0e-12 * max(xn, 1.0e-6) then
            break;
          end if;
        end for;
        xc := xn;
        zc := zn;
      end for;
      xH[i] := xc;
    end if;
  end for;
  // Helium (Saha, solved together with the electron density) and reionization.
  for i in 1:n loop
    xHe := 0.0;
    for it in 1:40 loop
      xHe := 0.5 * (xHe + heliumElectrons(z[i], (xH[i] + xHe) * nH0 * (1.0 + z[i])^3, par));
    end for;
    xe[i] := max(xH[i] + xHe, reionizationFraction(z[i], zReio, par));
    dtau[i] := xe[i] * nH0 * (1.0 + z[i])^3 * sigmaT * a[i] * GeneralRelativity.Constants.Mpc;
  end for;
  // Conformal time: eta(a1) by quadrature, then d(eta) = d(ln a) / (a H) along the grid.
  eta[1] := conformalTime(a[1], par);
  for i in 2:n loop
    eta[i] := eta[i - 1] + 0.5 * (a[i - 1] / a2H(a[i - 1], H0, Og + On, Ob + Oc, par.omegaK, OL) + a[i] / a2H(a[i], H0, Og + On, Ob + Oc, par.omegaK, OL))
      * log(a[i] / a[i - 1]);
  end for;
  tau[n] := 0.0;
  for i in n - 1:-1:1 loop
    tau[i] := tau[i + 1] + 0.5 * (dtau[i] + dtau[i + 1]) * (eta[i + 1] - eta[i]);
  end for;
  for i in 1:n loop
    g[i] := dtau[i] * exp(-tau[i]);
  end for;
  // Damping: k_D^-2 = int d(eta) / (6 (1 + R) dtau) [R^2 / (1 + R) + 16/15] (photon viscosity with polarization).
  sumD := 0.0;
  kD[1] := 1.0e30;
  for i in 2:n loop
    Ri := 0.75 * Ob / Og * a[i - 1];
    Rj := 0.75 * Ob / Og * a[i];
    wi := (Ri^2 / (1.0 + Ri) + 16.0 / 15.0) / (6.0 * (1.0 + Ri) * dtau[i - 1]);
    wj := (Rj^2 / (1.0 + Rj) + 16.0 / 15.0) / (6.0 * (1.0 + Rj) * dtau[i]);
    sumD := sumD + 0.5 * (wi + wj) * (eta[i] - eta[i - 1]);
    kD[i] := 1.0 / sqrt(max(sumD, 1.0e-300));
  end for;
  // Last scattering: the peak of g above z = 200 (the reionization bump is excluded), refined by a parabola.
  iMax := 2;
  for i in 2:n - 1 loop
    if z[i] > 200.0 and g[i] > g[iMax] then
      iMax := i;
    end if;
  end for;
  gl := g[iMax - 1];
  gm := g[iMax];
  gr := g[iMax + 1];
  shift := 0.5 * (gl - gr) / (gl - 2.0 * gm + gr);
  zStar := z[iMax] + shift * (if shift > 0.0 then z[iMax + 1] - z[iMax] else z[iMax] - z[iMax - 1]);
  etaStar := eta[iMax] + shift * (if shift > 0.0 then eta[iMax + 1] - eta[iMax] else eta[iMax] - eta[iMax - 1]);
  annotation(Documentation(info="<html>
<p>Hydrogen recombines in Saha equilibrium above <code>z = 1600</code> and then follows
<a href=\"modelica://GeneralRelativity.Cosmology.peeblesRate\">Peebles' equation</a>; helium is in Saha
equilibrium (it has recombined before hydrogen starts), and reionization adds CAMB's tanh profile at the
redshift that gives <code>par.tauReio</code>. From the free-electron density follow the scattering rate
<code>&tau;' = n<sub>e</sub>&sigma;<sub>T</sub>a</code>, the optical depth <code>&tau;</code> to today, the visibility
<code>g = &tau;'e<sup>&minus;&tau;</sup></code> (the probability density that a CMB photon last scattered at
conformal time <code>&eta;</code>; it integrates to one) and the diffusion-damping wavenumber
<code>k<sub>D</sub></code>, beyond which photon random walks have erased the acoustic waves.</p>
<p>With the Planck 2018 parameters: <code>z<sub>*</sub> = 1090</code>, <code>z<sub>re</sub> = 7.7</code>,
<code>k<sub>D</sub>(z<sub>*</sub>) &asymp; 0.14/Mpc</code>.</p>
</html>"));
end thermalHistory;
