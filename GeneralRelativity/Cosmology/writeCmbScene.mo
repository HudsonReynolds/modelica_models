within GeneralRelativity.Cosmology;
impure function writeCmbScene "Compute the CMB temperature spectrum and its ingredients, compare with Planck 2018, and write a viewer scene"
  input Parameters par;
  input String fileName "Scene file (full path)";
  input Integer nRec = 120 "Conformal-time samples across last scattering";
  input Integer nSnap = 40 "Snapshots of the transfer functions for the sound-wave picture";
  input Integer nEll = 80 "Logarithmically spaced multipoles from 30 to 2500 (2..29 are all computed)";
  output Real summary[24] "zStar, etaStar, eta0, rsStar, 100 thetaStar, kDStar, zReio, zDrag, l and D of peaks 1-3 (model), l and D of the first trough (model), chi2, nBins, rms vs best fit, peak 1 l and D (data), age";
protected
  constant Real pi = GeneralRelativity.Constants.pi;
  constant Integer nTh = 6000;
  constant Integer nLogK = 50;
  constant Integer nLinK = 300;
  constant Integer nHiK = 100;
  constant Integer nk = nLogK + nLinK + nHiK;
  constant Integer nIsw = 100;
  Real z[nTh], a[nTh], eta[nTh], xe[nTh], dtau[nTh], tau[nTh], g[nTh], kD[nTh];
  Real zStar, etaStar, zReio, eta0, rsStar, zDrag, tauDrag;
  Real H0, Og, On, Ob, Oc, OL;
  Real etaRec[nRec], etaIsw[nIsw], etaSnap[nSnap], aSnap[nSnap];
  Real gR[nRec], tauR[nRec], tauI[nIsw];
  constant Integer nTab = 3000;
  Real lnEtaStart, dLnEta, etaTab;
  Real dtauTab[nTab], kDTab[nTab], aTab[nTab];
  Real k[nk];
  Real rec[nk, nRec, 5], isw[nk, nIsw], snap[nk, nSnap, 5];
  Real ells[29 + nEll];
  Real Dl[29 + nEll], Dsw[29 + nEll], Ddop[29 + nEll], Disw[29 + nEll];
  Real M[29 + nEll], Msw[29 + nEll], Mdop[29 + nEll], Misw[29 + nEll];
  Real M2[29 + nEll, 4];
  Real lAll[2499], DAll[2499], swAll[2499], dopAll[2499], iswAll[2499];
  Real zLs[nTh], etaLs[nTh];
  Real eStart, eEnd;
  Real dataEll[111], dataD[111], dataLo[111], dataHi[111];
  Real bfEll[2507], bfD[2507];
  Real chi2, rms, nb, sig, modelValue;
  Real peaks[3, 2], trough[2], dataPeaks[3, 2], dataTrough[2];
  Integer np;
  Real lo, hi;
  Integer ib;
  Real yl, ym, yr, sh;
  Integer nRows;
  Real recRows[nTh, 6];
algorithm
  (H0, Og, On, Ob, Oc, OL) := densities(par);
  (z, a, eta, xe, dtau, tau, g, kD, zStar, etaStar, zReio) := thermalHistory(par, nTh);
  eta0 := conformalTime(1.0, par);
  rsStar := soundHorizon(1.0 / (1.0 + zStar), par);
  // Redshift-increasing copies for interpolation in z.
  for i in 1:nTh loop
    zLs[i] := z[nTh + 1 - i];
    etaLs[i] := eta[nTh + 1 - i];
  end for;
  // Baryon drag epoch: the optical depth of the baryons, int dtau / R from today, reaches 1.
  tauDrag := 0.0;
  zDrag := 0.0;
  for i in nTh - 1:-1:1 loop
    tauDrag := tauDrag + 0.5 * (dtau[i] / (0.75 * Ob / Og * a[i]) + dtau[i + 1] / (0.75 * Ob / Og * a[i + 1])) * (eta[i + 1] - eta[i]);
    if tauDrag >= 1.0 and zDrag == 0.0 then
      zDrag := z[i];
    end if;
  end for;
  // Conformal-time grids: across last scattering (z 1500 -> 700), then to today for the ISW term.
  eStart := Modelica.Math.Vectors.interpolate(zLs, etaLs, 1500.0);
  eEnd := Modelica.Math.Vectors.interpolate(zLs, etaLs, 700.0);
  for j in 1:nRec loop
    etaRec[j] := eStart + (eEnd - eStart) * (j - 1) / (nRec - 1);
    gR[j] := Modelica.Math.Vectors.interpolate(eta, g, etaRec[j]);
    tauR[j] := Modelica.Math.Vectors.interpolate(eta, tau, etaRec[j]);
  end for;
  for j in 1:30 loop
    etaIsw[j] := eEnd + (3.0 * etaStar - eEnd) * j / 30;
  end for;
  for j in 1:nIsw - 30 loop
    etaIsw[30 + j] := 3.0 * etaStar * ((eta0 - 1.0) / (3.0 * etaStar))^(j / (nIsw - 30));
  end for;
  for j in 1:nIsw loop
    tauI[j] := Modelica.Math.Vectors.interpolate(eta, tau, etaIsw[j]);
  end for;
  // Snapshots of the transfer functions, log-spaced in conformal time from a = 1e-6 to a = 0.02.
  for j in 1:nSnap loop
    etaSnap[j] := conformalTime(1.0e-6, par) * (conformalTime(0.02, par) / conformalTime(1.0e-6, par))^((j - 1) / (nSnap - 1));
    aSnap[j] := if etaSnap[j] < eta[1] then a2H(0.0, H0, Og + On, Ob + Oc, par.omegaK, OL) * etaSnap[j] + H0^2 * (Ob + Oc) / 4.0 * etaSnap[j]^2
      else Modelica.Math.Vectors.interpolate(eta, a, etaSnap[j]);
  end for;
  // Wavenumbers: logarithmic on the largest scales, uniform through the acoustic peaks, coarser beyond.
  for i in 1:nLogK loop
    k[i] := 1.0e-5 * (0.008 / 1.0e-5)^((i - 1) / nLogK);
  end for;
  for i in 1:nLinK loop
    k[nLogK + i] := 0.008 + (0.40 - 0.008) * (i - 1) / (nLinK - 1);
  end for;
  for i in 1:nHiK loop
    k[nLogK + nLinK + i] := 0.40 + 0.60 * i / nHiK;
  end for;
  // Thermal tables uniform in ln(eta) for the perturbation equations: scattering rate, damping scale, a.
  lnEtaStart := log(eta[1]);
  dLnEta := (log(eta0) - lnEtaStart) / (nTab - 1);
  for i in 1:nTab loop
    etaTab := exp(lnEtaStart + (i - 1) * dLnEta);
    dtauTab[i] := exp(Modelica.Math.Vectors.interpolate(eta, log(dtau), min(etaTab, eta[nTh])));
    kDTab[i] := Modelica.Math.Vectors.interpolate(eta, kD, min(etaTab, eta[nTh]));
    aTab[i] := Modelica.Math.Vectors.interpolate(eta, a, min(etaTab, eta[nTh]));
  end for;
  (rec, isw, snap) := acousticSources(par, k, etaRec, etaIsw, etaSnap, 3.0 * etaStar, eta0, lnEtaStart, dLnEta, dtauTab, kDTab, aTab,
    etaTightMax = Modelica.Math.Vectors.interpolate(zLs, etaLs, 2500.0));
  // Multipoles: every l below 30, then logarithmic to 2500.
  for i in 1:28 loop
    ells[i] := i + 1;
  end for;
  for i in 1:nEll + 1 loop
    ells[28 + i] := floor(30.0 * (2500.0 / 30.0)^((i - 1) / nEll) + 0.5);
  end for;
  (Dl, Dsw, Ddop, Disw) := cmbSpectrum(par, ells, k, etaRec, rec, etaIsw, isw, gR, tauR, tauI, eta0, dkFine = 1.0e-4);
  // Every integer l by cubic splines through the computed multipoles.
  M2 := splineSecondDerivatives(ells, transpose({Dl, Dsw, Ddop, Disw}));
  M := M2[:, 1];
  Msw := M2[:, 2];
  Mdop := M2[:, 3];
  Misw := M2[:, 4];
  for l in 2:2500 loop
    lAll[l - 1] := l;
    DAll[l - 1] := splineEvaluate(ells, Dl, M, l);
    swAll[l - 1] := max(splineEvaluate(ells, Dsw, Msw, l), 0.0);
    dopAll[l - 1] := max(splineEvaluate(ells, Ddop, Mdop, l), 0.0);
    iswAll[l - 1] := if l <= 400 then max(splineEvaluate(ells, Disw, Misw, l), 0.0) else 0.0;
  end for;
  // Peaks and the first trough of the model.
  np := 0;
  peaks := zeros(3, 2);
  for l in 100:1400 loop
    if np < 3 and DAll[l - 1] > DAll[l - 2] and DAll[l - 1] >= DAll[l] then
      np := np + 1;
      peaks[np, :] := {l, DAll[l - 1]};
    end if;
  end for;
  trough := {0.0, 1.0e30};
  for l in integer(peaks[1, 1]):integer(max(peaks[2, 1], peaks[1, 1] + 1)) loop
    if DAll[l - 1] < trough[2] then
      trough := {l, DAll[l - 1]};
    end if;
  end for;
  // Planck data: chi-square of the model (diagonal errors, model at the bin centres), the data's peaks.
  (dataEll, dataD, dataLo, dataHi) := Data.planck2018TT();
  (bfEll, bfD) := Data.planck2018BestFit();
  chi2 := 0.0;
  nb := 0;
  for i in 1:111 loop
    modelValue := splineEvaluate(ells, Dl, M, dataEll[i]);
    sig := if modelValue > dataD[i] then dataHi[i] else dataLo[i];
    chi2 := chi2 + ((modelValue - dataD[i]) / sig)^2;
    nb := nb + 1;
  end for;
  rms := 0.0;
  for l in 30:1500 loop
    rms := rms + (DAll[l - 1] / bfD[l - 1] - 1.0)^2;
  end for;
  rms := sqrt(rms / 1471.0);
  // Peaks of the binned data: highest bin in each window, refined by a parabola through its neighbours.
  for p in 1:3 loop
    lo := if p == 1 then 150.0 else if p == 2 then 450.0 else 720.0;
    hi := if p == 1 then 300.0 else if p == 2 then 620.0 else 900.0;
    ib := 0;
    for i in 29:110 loop
      if dataEll[i] >= lo and dataEll[i] <= hi and (ib == 0 or dataD[i] > dataD[max(ib, 1)]) then
        ib := i;
      end if;
    end for;
    yl := dataD[ib - 1];
    ym := dataD[ib];
    yr := dataD[ib + 1];
    sh := 0.5 * (yl - yr) / (yl - 2.0 * ym + yr);
    dataPeaks[p, :] := {dataEll[ib] + sh * 0.5 * (dataEll[ib + 1] - dataEll[ib - 1]), ym - 0.25 * (yl - yr) * sh};
  end for;
  summary := {zStar, etaStar, eta0, rsStar, 100.0 * rsStar / (eta0 - etaStar), Modelica.Math.Vectors.interpolate(eta, kD, etaStar), zReio, zDrag,
    peaks[1, 1], peaks[1, 2], peaks[2, 1], peaks[2, 2], peaks[3, 1], peaks[3, 2], trough[1], trough[2], chi2, nb, rms,
    dataPeaks[1, 1], dataPeaks[1, 2], cosmicTime(1.0, par), dataPeaks[2, 1], dataPeaks[3, 1]};
  // --- the scene ---------------------------------------------------------------------------------------
  GeneralRelativity.Visualization.beginScene(fileName, "CMB acoustic peaks (Planck 2018 LCDM)", "cartesian", 1.0, "Mpc", 1.0);
  GeneralRelativity.Visualization.writeView(fileName, "cmb");
  Modelica.Utilities.Streams.print("{\"type\":\"cosmology\",\"h\":" + String(par.h) + ",\"omegaB\":" + String(Ob) + ",\"omegaC\":" + String(Oc)
    + ",\"omegaM\":" + String(Ob + Oc) + ",\"omegaR\":" + String(Og + On) + ",\"omegaL\":" + String(OL) + ",\"omegaK\":" + String(par.omegaK)
    + ",\"Tcmb\":" + String(par.TCMB) + ",\"nS\":" + String(par.nS) + ",\"AS\":" + String(par.AS) + ",\"tauReio\":" + String(par.tauReio)
    + ",\"age\":" + String(summary[22]) + ",\"zEq\":" + String((Ob + Oc) / (Og + On) - 1.0) + ",\"zStar\":" + String(zStar)
    + ",\"zReio\":" + String(zReio) + ",\"zDrag\":" + String(zDrag) + ",\"etaStar\":" + String(etaStar) + ",\"eta0\":" + String(eta0)
    + ",\"rsStar\":" + String(rsStar) + ",\"thetaStar\":" + String(summary[5]) + ",\"kDStar\":" + String(summary[6])
    + ",\"DAstar\":" + String(eta0 - etaStar) + ",\"chi2\":" + String(chi2) + ",\"nBins\":" + String(nb) + "}", fileName);
  // Recombination history, about 300 rows from z = 3000 to today.
  nRows := 0;
  for i in 1:nTh loop
    if z[i] <= 3000.0 and mod(i, 10) == 0 then
      nRows := nRows + 1;
      recRows[nRows, :] := {z[i], xe[i], tau[i], g[i], eta[i], min(kD[i], 1.0e6)};
    end if;
  end for;
  Modelica.Utilities.Streams.print("{\"type\":\"recombination\",\"z\":" + GeneralRelativity.Visualization.jsonVector(recRows[1:nRows, 1], 6)
    + ",\"xe\":" + GeneralRelativity.Visualization.jsonVector(recRows[1:nRows, 2], 6)
    + ",\"tau\":" + GeneralRelativity.Visualization.jsonVector(recRows[1:nRows, 3], 6)
    + ",\"g\":" + GeneralRelativity.Visualization.jsonVector(recRows[1:nRows, 4], 6)
    + ",\"eta\":" + GeneralRelativity.Visualization.jsonVector(recRows[1:nRows, 5], 7)
    + ",\"kD\":" + GeneralRelativity.Visualization.jsonVector(recRows[1:nRows, 6], 6) + "}", fileName);
  Modelica.Utilities.Streams.print("{\"type\":\"cmbSpectrum\",\"ell\":" + GeneralRelativity.Visualization.jsonVector(lAll, 5)
    + ",\"Dl\":" + GeneralRelativity.Visualization.jsonVector(DAll, 6) + ",\"sw\":" + GeneralRelativity.Visualization.jsonVector(swAll, 6)
    + ",\"doppler\":" + GeneralRelativity.Visualization.jsonVector(dopAll, 6) + ",\"isw\":" + GeneralRelativity.Visualization.jsonVector(iswAll, 6)
    + ",\"computedEll\":" + GeneralRelativity.Visualization.jsonVector(ells, 5) + "}", fileName);
  Modelica.Utilities.Streams.print("{\"type\":\"cmbData\",\"source\":\"Planck 2018 (PR3) TT power spectrum: l < 30 unbinned (Commander), l >= 30 binned (Plik); Planck Collaboration, A&A 641, A1 and A6 (2020); Planck Legacy Archive via IRSA\""
    + ",\"ell\":" + GeneralRelativity.Visualization.jsonVector(dataEll, 7) + ",\"Dl\":" + GeneralRelativity.Visualization.jsonVector(dataD, 7)
    + ",\"errLo\":" + GeneralRelativity.Visualization.jsonVector(dataLo, 6) + ",\"errHi\":" + GeneralRelativity.Visualization.jsonVector(dataHi, 6) + "}", fileName);
  Modelica.Utilities.Streams.print("{\"type\":\"cmbTheory\",\"source\":\"Planck 2018 best-fit base LCDM (TT,TE,EE+lowE+lensing), lensed, CAMB\""
    + ",\"ell\":" + GeneralRelativity.Visualization.jsonVector(bfEll, 5) + ",\"Dl\":" + GeneralRelativity.Visualization.jsonVector(bfD, 6) + "}", fileName);
  Modelica.Utilities.Streams.print("{\"type\":\"cmbComparison\",\"chi2\":" + String(chi2) + ",\"nBins\":" + String(nb) + ",\"rmsVsTheory\":" + String(rms)
    + ",\"peaksModel\":" + GeneralRelativity.Visualization.jsonMatrix(peaks, 6) + ",\"troughModel\":" + GeneralRelativity.Visualization.jsonVector(trough, 6)
    + ",\"peaksData\":" + GeneralRelativity.Visualization.jsonMatrix(dataPeaks, 6) + "}", fileName);
  Modelica.Utilities.Streams.print("{\"type\":\"acoustic\",\"k\":" + GeneralRelativity.Visualization.jsonVector(k, 6)
    + ",\"eta\":" + GeneralRelativity.Visualization.jsonVector(etaSnap, 6) + ",\"a\":" + GeneralRelativity.Visualization.jsonVector(aSnap, 6)
    + ",\"z\":" + GeneralRelativity.Visualization.jsonVector({1.0 / aSnap[j] - 1.0 for j in 1:nSnap}, 6)
    + ",\"deltaC\":" + GeneralRelativity.Visualization.jsonMatrix(transpose(snap[:, :, 1]), 5)
    + ",\"deltaB\":" + GeneralRelativity.Visualization.jsonMatrix(transpose(snap[:, :, 2]), 5)
    + ",\"deltaG\":" + GeneralRelativity.Visualization.jsonMatrix(transpose(snap[:, :, 3]), 5)
    + ",\"deltaNu\":" + GeneralRelativity.Visualization.jsonMatrix(transpose(snap[:, :, 4]), 5)
    + ",\"phi\":" + GeneralRelativity.Visualization.jsonMatrix(transpose(snap[:, :, 5]), 5) + "}", fileName);
  annotation(Documentation(info="<html>
<p>The work of <a href=\"modelica://GeneralRelativity.Examples.CosmicMicrowaveBackground\">Examples.CosmicMicrowaveBackground</a>.
Records written (one JSON object per line):</p>
<ul>
<li><code>cosmology</code>: parameters and derived quantities (<code>zStar, etaStar, eta0, rsStar, thetaStar</code> (100&theta;<sub>*</sub>),
<code>kDStar, DAstar, zDrag, zReio, chi2, nBins</code>, ...); conformal quantities in Mpc;</li>
<li><code>recombination</code>: <code>z, xe, tau, g, eta, kD</code> from <code>z = 3000</code> to 0;</li>
<li><code>cmbSpectrum</code>: <code>D<sub>l</sub></code> (&mu;K&sup2;) for every <code>l = 2..2500</code> (spline through the computed
<code>computedEll</code>) and the single-term spectra <code>sw, doppler, isw</code>;</li>
<li><code>cmbData</code>, <code>cmbTheory</code>: Planck 2018 measurements and best-fit spectrum;</li>
<li><code>cmbComparison</code>: chi-square, rms difference from the best fit, peak positions and heights;</li>
<li><code>acoustic</code>: transfer functions <code>deltaC, deltaB, deltaG, deltaNu, phi</code> per unit primordial
curvature, rows = the snapshot times <code>eta, a, z</code>, columns = the wavenumbers <code>k</code> (1/Mpc,
<code>10<sup>&minus;5</sup>..1</code>; uniform at 0.0013/Mpc from 0.008 to 0.4, 0.006/Mpc beyond).</li>
</ul>
</html>"));
end writeCmbScene;
