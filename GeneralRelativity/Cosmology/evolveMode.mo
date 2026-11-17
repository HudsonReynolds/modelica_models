within GeneralRelativity.Cosmology;
function evolveMode "Evolve one Fourier mode from adiabatic initial conditions (unit primordial curvature) and sample it"
  input Real k(unit = "1/Mpc");
  input Real etaRec[:] "Increasing conformal times for the recombination sources";
  input Real etaIsw[:] "Increasing conformal times for the integrated Sachs-Wolfe source";
  input Real etaSnap[:] "Increasing conformal times for snapshots of the transfer functions";
  input Real etaOff "Radiation perturbations of deep sub-horizon modes are dropped after max(etaOff, 50/k)";
  input Real eta0 "Today";
  input Real lnEtaStart "Thermal tables: ln(eta) of the first entry";
  input Real dLnEta "Thermal tables: spacing in ln(eta)";
  input Real dtauTable[:] "Thomson scattering rate n_e sigma_T a, 1/Mpc";
  input Real kDTable[size(dtauTable, 1)] "Diffusion-damping wavenumber, 1/Mpc";
  input Real aTable[size(dtauTable, 1)] "Scale factor";
  input Real H0;
  input Real omegaGamma;
  input Real omegaNu;
  input Real omegaB;
  input Real omegaC;
  input Real omegaK;
  input Real omegaLambda;
  input Real rtol = 1.0e-6;
  input Real atol = 1.0e-10;
  input Real tightCoupling = 30.0 "The photon hierarchy starts where the scattering rate falls below this times max(k, a'/a)";
  input Real etaTightMax = 1.0e30 "... and at the latest here (before recombination, so no mode changes equations during last scattering)";
  output Real rec[size(etaRec, 1), 5] "delta_gamma/4, psi, theta_b/k, psi' + phi', Pi = F_g2 + G_g0 + G_g2";
  output Real isw[size(etaIsw, 1)] "psi' + phi'";
  output Real snap[size(etaSnap, 1), 5] "delta_c, delta_b, delta_gamma, delta_nu, phi";
  output Real etaTight "End of tight coupling for this mode";
protected
  // Dormand-Prince 5(4) coefficients.
  constant Real c2 = 1.0 / 5.0, c3 = 3.0 / 10.0, c4 = 4.0 / 5.0, c5 = 8.0 / 9.0;
  constant Real a21 = 1.0 / 5.0;
  constant Real a31 = 3.0 / 40.0, a32 = 9.0 / 40.0;
  constant Real a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0;
  constant Real a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0, a53 = 64448.0 / 6561.0, a54 = -212.0 / 729.0;
  constant Real a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0, a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0, a65 = -5103.0 / 18656.0;
  constant Real b1 = 35.0 / 384.0, b3 = 500.0 / 1113.0, b4 = 125.0 / 192.0, b5 = -2187.0 / 6784.0, b6 = 11.0 / 84.0;
  constant Real e1 = 71.0 / 57600.0, e3 = -71.0 / 16695.0, e4 = 71.0 / 1920.0, e5 = -17253.0 / 339200.0, e6 = 22.0 / 525.0, e7 = -1.0 / 40.0;
  constant Integer n = 35;
  Integer nr = size(etaRec, 1);
  Integer ni = size(etaIsw, 1);
  Integer ns = size(etaSnap, 1);
  Integer nt = size(dtauTable, 1);
  Integer ir = 1;
  Integer ii = 1;
  Integer is = 1;
  Real y[n];
  Real yn[n];
  Real k1[n], k2[n], k3[n], k4[n], k5[n], k6[n], k7[n];
  Real err[n];
  Real d[n];
  Real eta;
  Real etaInit;
  Real etaNext;
  Real etaStop;
  Real h;
  Real errNorm;
  Real sc;
  Integer phase;
  Real Rnu = omegaNu / (omegaGamma + omegaNu);
  Real psi0;
  Real a;
  Real psiOf;
  Real dIsw = 0.0;
  Real eOff = min(max(etaOff, 50.0 / k), eta0);
  Real etaT;
  Real aT;
  Real D;
  Real R;
  Real theta;
  Real dgOld;
  Boolean accept;
algorithm
  rec := zeros(nr, 5);
  isw := zeros(ni);
  snap := zeros(ns, 5);
  // End of tight coupling: the first table entry where dtau < tightCoupling * max(k, a'/a).
  etaTight := min(eOff, etaTightMax);
  for i in 1:nt loop
    etaT := exp(lnEtaStart + (i - 1) * dLnEta);
    aT := aTable[i];
    if dtauTable[i] < tightCoupling * max(k, a2H(aT, H0, omegaGamma + omegaNu, omegaB + omegaC, omegaK, omegaLambda) / aT) then
      etaTight := min(etaT, etaTight);
      break;
    end if;
  end for;
  // Adiabatic initial conditions deep in the radiation era (Ma & Bertschinger eq. 96 with C = 1/2,
  // i.e. unit comoving curvature): a from the exact radiation + matter solution.
  etaInit := min(1.0e-3 / k, 0.01);
  psi0 := 10.0 / (15.0 + 4.0 * Rnu);
  y := zeros(n);
  y[1] := H0 * sqrt(omegaGamma + omegaNu) * etaInit + H0^2 * (omegaB + omegaC) / 4.0 * etaInit^2;
  y[2] := (1.0 + 0.4 * Rnu) * psi0;
  y[3] := -1.5 * psi0;
  y[5] := -1.5 * psi0;
  y[4] := 0.5 * k^2 * etaInit * psi0;
  y[6] := y[4];
  y[7] := -2.0 * psi0;
  y[8] := y[4];
  y[9] := -2.0 * psi0;
  y[10] := y[4];
  y[11] := 2.0 * (k * etaInit)^2 * psi0 / 15.0;
  eta := etaInit;
  h := 0.1 * etaInit;
  phase := 0;
  while eta < eta0 and (ir <= nr or ii <= ni or is <= ns) loop
    // Next stop: an output time or a phase boundary.
    etaNext := eta0;
    if ir <= nr then
      etaNext := min(etaNext, etaRec[ir]);
    end if;
    if ii <= ni then
      etaNext := min(etaNext, etaIsw[ii]);
    end if;
    if is <= ns then
      etaNext := min(etaNext, etaSnap[is]);
    end if;
    etaStop := if phase == 0 and etaTight < etaNext then etaTight else if phase == 1 and eOff < etaNext then eOff else etaNext;
    // Adaptive Dormand-Prince steps up to etaStop.
    k1 := modeDerivatives(y, eta, k, phase, tableValue(lnEtaStart, dLnEta, dtauTable, eta), H0, omegaGamma, omegaNu, omegaB, omegaC, omegaK, omegaLambda);
    while eta < etaStop loop
      h := min(h, etaStop - eta);
      k2 := modeDerivatives(y + h * a21 * k1, eta + c2 * h, k, phase, tableValue(lnEtaStart, dLnEta, dtauTable, eta + c2 * h), H0, omegaGamma, omegaNu, omegaB, omegaC, omegaK, omegaLambda);
      k3 := modeDerivatives(y + h * (a31 * k1 + a32 * k2), eta + c3 * h, k, phase, tableValue(lnEtaStart, dLnEta, dtauTable, eta + c3 * h), H0, omegaGamma, omegaNu, omegaB, omegaC, omegaK, omegaLambda);
      k4 := modeDerivatives(y + h * (a41 * k1 + a42 * k2 + a43 * k3), eta + c4 * h, k, phase, tableValue(lnEtaStart, dLnEta, dtauTable, eta + c4 * h), H0, omegaGamma, omegaNu, omegaB, omegaC, omegaK, omegaLambda);
      k5 := modeDerivatives(y + h * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4), eta + c5 * h, k, phase, tableValue(lnEtaStart, dLnEta, dtauTable, eta + c5 * h), H0, omegaGamma, omegaNu, omegaB, omegaC, omegaK, omegaLambda);
      k6 := modeDerivatives(y + h * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5), eta + h, k, phase, tableValue(lnEtaStart, dLnEta, dtauTable, eta + h), H0, omegaGamma, omegaNu, omegaB, omegaC, omegaK, omegaLambda);
      yn := y + h * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6);
      k7 := modeDerivatives(yn, eta + h, k, phase, tableValue(lnEtaStart, dLnEta, dtauTable, eta + h), H0, omegaGamma, omegaNu, omegaB, omegaC, omegaK, omegaLambda);
      err := h * (e1 * k1 + e3 * k3 + e4 * k4 + e5 * k5 + e6 * k6 + e7 * k7);
      errNorm := 0.0;
      for i in 1:n loop
        sc := atol + rtol * max(abs(y[i]), abs(yn[i]));
        errNorm := errNorm + (err[i] / sc)^2;
      end for;
      errNorm := sqrt(errNorm / n);
      accept := errNorm <= 1.0;
      if accept then
        eta := if etaStop - (eta + h) < 1.0e-12 * etaStop then etaStop else eta + h;
        y := yn;
        k1 := k7;
      end if;
      h := h * min(5.0, max(0.2, 0.9 * (max(errNorm, 1.0e-10))^(-0.2)));
    end while;
    // Phase changes.
    if phase == 0 and eta >= etaTight then
      // Leave tight coupling: the photon shear and polarization take their tight-coupling values
      // (sigma = 16/45 theta/dtau, G_0 = 5/4 F_2, G_2 = 1/4 F_2), and the diffusion damping accumulated so
      // far, exp(-(k/k_D)^2), is applied to the oscillating part Theta_0 + (1 + R) psi and to the velocity.
      a := y[1];
      R := 0.75 * omegaB / omegaGamma * a;
      D := exp(-(k / tableValue(lnEtaStart, dLnEta, kDTable, eta))^2);
      psiOf := y[2] - 3.0 * H0^2 * omegaNu * y[11] / (a^2 * k^2);
      dgOld := y[7];
      y[7] := D * (y[7] + 4.0 * (1.0 + R) * psiOf) - 4.0 * (1.0 + R) * psiOf;
      y[5] := y[5] + 0.75 * (y[7] - dgOld);
      theta := D * y[8];
      y[6] := theta;
      y[8] := theta;
      y[16] := 32.0 / 45.0 * theta / tableValue(lnEtaStart, dLnEta, dtauTable, eta);
      y[25] := 1.25 * y[16];
      y[27] := 0.25 * y[16];
      phase := 1;
      h := min(h, 0.1 / tableValue(lnEtaStart, dLnEta, dtauTable, eta));
    end if;
    if phase == 1 and eta >= eOff then
      phase := 2;
    end if;
    // Outputs due now.
    a := y[1];
    psiOf := y[2] - (if phase < 2 then 3.0 * H0^2 * (omegaNu * y[11] + (if phase == 1 then omegaGamma * y[16] else 0.0)) / (a^2 * k^2) else 0.0);
    if (ir <= nr and etaRec[ir] <= eta * (1.0 + 1.0e-12)) or (ii <= ni and etaIsw[ii] <= eta * (1.0 + 1.0e-12)) then
      // psi' + phi', the integrated Sachs-Wolfe source (psi' from the anisotropic-stress term).
      d := modeDerivatives(y, eta, k, phase, tableValue(lnEtaStart, dLnEta, dtauTable, eta), H0, omegaGamma, omegaNu, omegaB, omegaC, omegaK, omegaLambda);
      dIsw := d[2] + (if phase < 2 then d[2] - 3.0 * H0^2 / k^2 * ((omegaNu * d[11] + (if phase == 1 then omegaGamma * d[16] else 0.0)) / a^2
        - 2.0 * (omegaNu * y[11] + (if phase == 1 then omegaGamma * y[16] else 0.0)) * d[1] / a^3) else d[2]);
    end if;
    while ir <= nr and etaRec[ir] <= eta * (1.0 + 1.0e-12) loop
      rec[ir, :] := {y[7] / 4.0, psiOf, y[6] / k, dIsw, if phase == 1 then y[16] + y[25] + y[27] else 0.0};
      ir := ir + 1;
    end while;
    while ii <= ni and etaIsw[ii] <= eta * (1.0 + 1.0e-12) loop
      isw[ii] := dIsw;
      ii := ii + 1;
    end while;
    while is <= ns and etaSnap[is] <= eta * (1.0 + 1.0e-12) loop
      snap[is, :] := {y[3], y[5], y[7], y[9], y[2]};
      is := is + 1;
    end while;
  end while;
  annotation(Documentation(info="<html>
<p>Adaptive Dormand-Prince 5(4) integration of
<a href=\"modelica://GeneralRelativity.Cosmology.modeDerivatives\">modeDerivatives</a> in conformal time from
<code>k&eta; = 10<sup>&minus;3</sup></code>, stepping exactly onto every requested output time and phase
boundary. Initial conditions are adiabatic and normalized to unit comoving curvature
<code>&Rscr;</code>, so outputs are transfer functions: multiply by the primordial amplitude.</p>
<p>Photons and baryons start tightly coupled; the photon hierarchy takes over where the scattering rate
<code>&tau;'</code> falls below <code>tightCoupling &times; max(k, a'/a)</code> (early for small scales, during
recombination for large ones), and at the latest at <code>etaTightMax</code>, chosen before recombination so that
no mode changes equations while the CMB is emitted. The stiffness of the hierarchy (rate
<code>&tau;'(1 + 1/R)</code>) sets the step size there.</p>
</html>"));
end evolveMode;
