within GeneralRelativity.Visualization;
impure function writeLensInversion "Invert a lensed image (Lensing.Inversion.invert) and append the case to a scene for the viewer's lensfit style: observed, first and fitted simulations, residual and fit summary"
  input String fileName "Scene file, already started with beginScene";
  input String caseName = "" "Case key when one scene holds several inversions; empty for a single one";
  input String title = "" "Title shown by the viewer";
  input String description = "" "One or two sentences on what the image shows";
  input Real u[:](each unit = "rad") "Pixel-centre angles along the first axis";
  input Real v[:](each unit = "rad") "Pixel-centre angles along the second axis";
  input Real observed[size(u, 1), size(v, 1)] "The image, lens light subtracted";
  input Real noiseSigma "Pixel noise standard deviation";
  input Boolean pointMass = false;
  input Real psfSigma(unit = "rad") = 0.0 "Gaussian seeing";
  input Real distances[3](each unit = "m") "D_l, D_s, D_ls";
  input Real apertureM(unit = "m") = 10.0e-3 * GeneralRelativity.Constants.Mpc "Aperture for the galaxy mass reported to the lensing style";
  input Real trueLens[:] = fill(0.0, 0) "The true lens (7 entries) when the image is synthetic; empty when unknown";
  input Real trueSource[:] = fill(0.0, 0) "The true source (7 entries); empty when unknown";
  input Integer maxIterations = 30;
  input Integer startIterations = 4 "Iterations of each trial start in invert";
  input Real threshold = 0.2;
  input Boolean writeLensRecord = false "Also write a lens record of the fit, which the lensing style draws a sky through (one per scene)";
  input Real zLens = 0.5 "For the lens record";
  input Real zSource = 2.0 "For the lens record";
  input Real H0 = 70.0 "For the lens record";
  input Real OmegaM = 0.3 "For the lens record";
  output Real ringRadius(unit = "rad") "Einstein radius read off the image (step 2)";
  output Real massInitialSolar "Mass inside that radius (step 2)";
  output Real lens[7] "Fitted lens";
  output Real sources[1, 7] "Fitted source";
  output Real lensSigma[7];
  output Real thetaE(unit = "rad") "Fitted equal-area Einstein radius";
  output Real thetaESigma(unit = "rad");
  output Real massSolar "Fitted projected mass inside the critical curve";
  output Real massSigmaSolar;
  output Real massInApertureSolar "Mass within the aperture implied by the fitted profile";
  output Real trueThetaE(unit = "rad") "Equal-area Einstein radius of the true lens; 0 when unknown";
  output Real trueMassSolar "Mass inside the true critical curve; 0 when unknown";
  output Real chi2Initial "Cost of the first simulation (step 4)";
  output Real chi2Final "Cost of the fitted simulation";
  output Real explained;
  output Real correlation;
  output Integer iterations;
  output Boolean converged;
protected
  constant Real arcsec = GeneralRelativity.Constants.arcsec;
  constant Real deg = GeneralRelativity.Constants.pi / 180.0;
  constant Real GMsun = GeneralRelativity.Constants.GM_sun;
  Integer nHistory = startIterations + maxIterations + 1;
  Real ringCentre[2];
  Real lens0[7];
  Real sources0[1, 7];
  Real chi2History[startIterations + maxIterations + 1];
  Real lensHistory[startIterations + maxIterations + 1, 7];
  Real massHistorySolar[startIterations + maxIterations + 1];
  Real sourcesSigma[1, 7];
  Real lensCovariance[7, 7];
  Real firstGuess[size(u, 1), size(v, 1)];
  Real fitted[size(u, 1), size(v, 1)];
  Real sumSquares;
  Real GM;
  Real GMSigma;
  Real thetaHistory;
  Real sigmaHistory;
  Real GMHistory;
  Real GMHistorySigma;
  Real trueSigma;
  Real trueGM;
  Real trueGMSigma;
  Real toDisplay[14] = {1.0 / arcsec, 1.0 / arcsec, 1.0 / arcsec, 1.0, 1.0 / deg, 1.0, 1.0 / deg, 1.0 / arcsec, 1.0 / arcsec, 1.0, 1.0 / arcsec, 1.0, 1.0 / deg, 1.0};
  Boolean known = size(trueLens, 1) == 7 and size(trueSource, 1) == 7;
  Real truthDisplay[14] "The truth in display units; zeros when unknown";
algorithm
  (ringRadius, ringCentre, lens0, sources0, lens, sources, chi2History, lensHistory, iterations, converged, lensSigma, sourcesSigma, lensCovariance) :=
    GeneralRelativity.Lensing.Inversion.invert(u, v, observed, noiseSigma, pointMass, psfSigma, maxIterations, startIterations, threshold);
  massInitialSolar := GeneralRelativity.Lensing.massInEinsteinRadius(ringRadius, distances[1], distances[2], distances[3]) / GMsun;
  firstGuess := GeneralRelativity.Lensing.Inversion.renderImage(u, v, lens0, sources0, pointMass, 2, psfSigma);
  (sumSquares, chi2Initial, explained, correlation) := GeneralRelativity.Lensing.Inversion.imageResidual(observed, firstGuess, noiseSigma);
  fitted := GeneralRelativity.Lensing.Inversion.renderImage(u, v, lens, sources, pointMass, 2, psfSigma);
  (sumSquares, chi2Final, explained, correlation) := GeneralRelativity.Lensing.Inversion.imageResidual(observed, fitted, noiseSigma);
  (thetaE, thetaESigma, GM, GMSigma) := GeneralRelativity.Lensing.Inversion.einsteinMass(lens, lensCovariance, pointMass, distances[1], distances[2], distances[3]);
  massSolar := GM / GMsun;
  massSigmaSolar := GMSigma / GMsun;
  massInApertureSolar := if pointMass then massSolar else massSolar * apertureM / (distances[1] * thetaE);
  for k in 1:nHistory loop
    (thetaHistory, sigmaHistory, GMHistory, GMHistorySigma) := GeneralRelativity.Lensing.Inversion.einsteinMass(lensHistory[k, :], zeros(7, 7), pointMass, distances[1], distances[2], distances[3]);
    massHistorySolar[k] := GMHistory / GMsun;
  end for;
  if known then
    (trueThetaE, trueSigma, trueGM, trueGMSigma) := GeneralRelativity.Lensing.Inversion.einsteinMass(trueLens, zeros(7, 7), pointMass, distances[1], distances[2], distances[3]);
    trueMassSolar := trueGM / GMsun;
    truthDisplay := cat(1, trueLens[1:7], trueSource[1:7]) .* toDisplay;
  else
    trueThetaE := 0.0;
    trueMassSolar := 0.0;
    truthDisplay := zeros(14);
  end if;
  writeLensImage(fileName, "observed", u / arcsec, v / arcsec, observed, "", caseName);
  writeLensImage(fileName, "initial", u / arcsec, v / arcsec, firstGuess, "", caseName);
  writeLensImage(fileName, "model", u / arcsec, v / arcsec, fitted, "", caseName);
  writeLensImage(fileName, "residual", u / arcsec, v / arcsec, (observed - fitted) / noiseSigma, "sigma", caseName);
  writeLensFit(fileName, if pointMass then "point" else "galaxy",
    {"b, deflection scale (arcsec)", "lens x (arcsec)", "lens y (arcsec)", "axis ratio q", "position angle (deg)", "shear", "shear angle (deg)",
     "source x (arcsec)", "source y (arcsec)", "source I_e", "source r_e (arcsec)", "source q", "source angle (deg)", "Sersic n"},
    cat(1, lens0, sources0[1, :]) .* toDisplay, cat(1, lens, sources[1, :]) .* toDisplay, cat(1, lensSigma, sourcesSigma[1, :]) .* toDisplay,
    truthDisplay[1:(if known then 14 else 0)], chi2History, massHistorySolar, massSolar, massSigmaSolar,
    massInitialSolar, trueMassSolar, thetaE / arcsec, distances / GeneralRelativity.Constants.Mpc, iterations, converged, size(u, 1) * size(v, 1),
    explained, correlation, caseName, title, description);
  if writeLensRecord then
    // The lensing style sets b from a mass within the aperture through the isothermal-sphere relation, which is this.
    writeLens(fileName, if pointMass then "point" else "galaxy",
      if pointMass then massSolar else GeneralRelativity.Lensing.massInEinsteinRadius(lens[1], distances[1], distances[2], distances[3]) / GMsun * apertureM / (distances[1] * lens[1]),
      apertureM / (1.0e-3 * GeneralRelativity.Constants.Mpc), zLens, zSource, H0, OmegaM, lens[4], lens[5] / deg, lens[6], lens[7] / deg, lens[1] / arcsec,
      distances / GeneralRelativity.Constants.Mpc);
  end if;
  annotation(Documentation(info="<html>
<p>Runs <a href=\"modelica://GeneralRelativity.Lensing.Inversion.invert\">Lensing.Inversion.invert</a> on one image, measures the mass
with <a href=\"modelica://GeneralRelativity.Lensing.Inversion.einsteinMass\">einsteinMass</a>, and appends four
<code>lensImage</code> records (observed, first simulation, fitted simulation, residual in units of the noise) and a
<code>lensFit</code> record, all tagged with <code>caseName</code>. Call it once per image to build a library that the viewer's
<b>lensfit</b> style shows as a gallery
(<a href=\"modelica://GeneralRelativity.Examples.LensInversionLibrary\">Examples.LensInversionLibrary</a>). When the true lens is
given, the fit record carries it and the mass inside the true critical curve, and the viewer grades the fit against them.</p>
</html>"));
end writeLensInversion;
