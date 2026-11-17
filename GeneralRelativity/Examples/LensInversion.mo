within GeneralRelativity.Examples;
model LensInversion "Lens inversion on a synthetic image: read the Einstein radius off the arcs, estimate the mass, simulate, compare and iterate until the simulation matches; the true mass is known and the recovery is asserted"
  // The truth, used only to make the image and to grade the answer.
  parameter Real trueMassSolar = 5.0e11 "True lens mass in solar masses: projected within the aperture (galaxy) or total (point)";
  parameter Real apertureKpc = 10.0 "Aperture radius that defines the galaxy mass";
  parameter Real zLens = 0.5 "Lens redshift";
  parameter Real zSource = 2.0 "Source redshift";
  parameter Real H0 = 70.0 "Hubble constant, km/s/Mpc";
  parameter Real OmegaM = 0.3 "Matter density (flat universe)";
  parameter Boolean pointMass = false "true: point-mass lens; false: isothermal ellipse with external shear";
  parameter Real axisRatio = 0.75 "True galaxy axis ratio";
  parameter Real positionAngleDeg = 30.0 "True galaxy position angle";
  parameter Real shear = 0.05 "True external shear";
  parameter Real shearAngleDeg = 110.0 "True shear position angle";
  parameter Real lensCentreArcsec[2] = {0.10, -0.05} "True lens position in the image";
  parameter Real sourceArcsec[2] = {0.22, -0.14} "True source position behind the lens";
  parameter Real sourceRadiusArcsec = 0.18 "True source half-light radius";
  parameter Real sourceAxisRatio = 0.7;
  parameter Real sourceAngleDeg = 50.0;
  parameter Real sersicIndex = 1.0 "True Sersic index of the source";
  // The observation.
  parameter Integer pixels = 64 "Image size per axis";
  parameter Real fieldArcsec = 3.5 "Half-width of the image";
  parameter Real psfArcsec = 0.06 "Gaussian point-spread width (seeing), assumed known";
  parameter Real noiseFraction = 0.02 "Pixel noise as a fraction of the brightest noiseless pixel";
  parameter Integer seed = 7 "Noise seed";
  // The inversion.
  parameter Real threshold = 0.2 "Fraction of the peak below which pixels do not count for the centroid and the ring radius";
  parameter Integer maxIterations = 30;
  parameter Boolean checkRecovery = true "false switches off the assertions that grade the fit against the truth";
  parameter String sceneFile = "lens_inversion.ndjson" "Scene file name";
  parameter String sceneDirectory = "" "Output folder; empty uses GR_SCENE_DIR, else ~/Documents/Visualization";
  parameter String imageFile = "lens_inversion_observed.mat" "The observed image is also saved here as the bordered matrix lensImage (MATLAB v4), the input of LensInversionFromFile; empty to skip";
  final parameter Real arcsec = GeneralRelativity.Constants.arcsec;
  final parameter Real deg = GeneralRelativity.Constants.pi / 180.0;
  final parameter Real Dl(unit = "m") = GeneralRelativity.Lensing.angularDiameterDistance(0.0, zLens, H0, OmegaM);
  final parameter Real Ds(unit = "m") = GeneralRelativity.Lensing.angularDiameterDistance(0.0, zSource, H0, OmegaM);
  final parameter Real Dls(unit = "m") = GeneralRelativity.Lensing.angularDiameterDistance(zLens, zSource, H0, OmegaM);
  final parameter Real apertureM(unit = "m") = apertureKpc * 1.0e-3 * GeneralRelativity.Constants.Mpc;
  final parameter Real trueGM(unit = "m3/s2") = trueMassSolar * GeneralRelativity.Constants.GM_sun;
  final parameter Real trueThetaE(unit = "rad") = if pointMass then GeneralRelativity.Lensing.einsteinRadiusPointMass(trueGM, Dl, Ds, Dls)
    else GeneralRelativity.Lensing.einsteinRadiusSIS(sqrt(trueGM / (GeneralRelativity.Constants.pi * apertureM)), Ds, Dls) "True Einstein radius";
  final parameter Real trueMassInRingSolar = GeneralRelativity.Lensing.massInEinsteinRadius(trueThetaE, Dl, Ds, Dls) / GeneralRelativity.Constants.GM_sun
    "True projected mass inside the Einstein radius: the whole mass of a point, theta_E D_l / aperture of the galaxy mass";
  final parameter Real trueLens[7] = {trueThetaE, lensCentreArcsec[1] * arcsec, lensCentreArcsec[2] * arcsec, axisRatio, positionAngleDeg * deg, shear, shearAngleDeg * deg};
  final parameter Real trueSources[1, 7] = [sourceArcsec[1] * arcsec, sourceArcsec[2] * arcsec, 1.0, sourceRadiusArcsec * arcsec, sourceAxisRatio, sourceAngleDeg * deg, sersicIndex];
  final parameter Real u[pixels](each unit = "rad") = {(-fieldArcsec + (i - 0.5) * 2.0 * fieldArcsec / pixels) * arcsec for i in 1:pixels} "Pixel centres";
  final parameter Real v[pixels](each unit = "rad") = u;
  final parameter Real toDisplay[14] = {1.0 / arcsec, 1.0 / arcsec, 1.0 / arcsec, 1.0, 1.0 / deg, 1.0, 1.0 / deg, 1.0 / arcsec, 1.0 / arcsec, 1.0, 1.0 / arcsec, 1.0, 1.0 / deg, 1.0}
    "Radians to arcseconds and degrees for the scene";
  // Step 1: the image.
  parameter Real noiseless[pixels, pixels](each fixed = false);
  parameter Real noiseSigma(fixed = false) "Pixel noise standard deviation";
  parameter Real observed[pixels, pixels](each fixed = false) "The input image";
  // Step 2: the Einstein radius read off the image and the mass it implies.
  parameter Real centre[2](each fixed = false, each unit = "rad") "Lens position from the circle fitted to the bright pixels";
  parameter Real thetaE0(fixed = false, unit = "rad") "Einstein radius of the circle fitted to the bright pixels";
  parameter Real ringSpread(fixed = false, unit = "rad");
  parameter Real massInitialSolar(fixed = false) "Projected mass inside thetaE0, in solar masses: the radius-based estimate";
  // Step 3: the simulation of that estimate.
  parameter Real lens0[7](each fixed = false);
  parameter Real sources0[1, 7](each fixed = false);
  parameter Real firstGuess[pixels, pixels](each fixed = false) "Simulated image of the first estimate";
  // Step 4: comparison.
  parameter Real chi2Initial(fixed = false);
  parameter Real explainedInitial(fixed = false);
  parameter Real correlationInitial(fixed = false);
  // Step 5: iterate until the simulation matches.
  parameter Real lens[7](each fixed = false) "Fitted lens";
  parameter Real sources[1, 7](each fixed = false) "Fitted source";
  parameter Real chi2History[maxIterations + 1](each fixed = false);
  parameter Real lensHistory[maxIterations + 1, 7](each fixed = false);
  parameter Integer iterations(fixed = false);
  parameter Boolean converged(fixed = false);
  parameter Real lensSigma[7](each fixed = false) "Standard errors of the lens parameters";
  parameter Real sourcesSigma[1, 7](each fixed = false);
  parameter Real fitted[pixels, pixels](each fixed = false) "Simulated image of the fit";
  parameter Real chi2Final(fixed = false) "Should be about pixels^2 when the model and the noise estimate are right";
  parameter Real explainedFinal(fixed = false);
  parameter Real correlationFinal(fixed = false);
  parameter Real massSolar(fixed = false) "Fitted projected mass inside the Einstein radius, in solar masses";
  parameter Real massSigmaSolar(fixed = false) "Its standard error (the mass goes as thetaE^2)";
  parameter Real massInApertureSolar(fixed = false) "The galaxy's mass within the aperture implied by the fit, comparable to trueMassSolar";
  parameter Real massHistorySolar[maxIterations + 1](each fixed = false);
  parameter Real massErrorSigmas(fixed = false) "(fitted - true mass inside the ring) / standard error";
protected
  parameter Real residual[pixels, pixels](each fixed = false);
  parameter Real sumSquaresInitial(fixed = false);
  parameter Real sumSquaresFinal(fixed = false);
  parameter Real initialDisplay[14](each fixed = false);
  parameter Real fittedDisplay[14](each fixed = false);
  parameter Real sigmaDisplay[14](each fixed = false);
  parameter Real truthDisplay[14](each fixed = false);
  parameter String path = GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory);
  parameter String imagePath = if imageFile == "" then "" else GeneralRelativity.Visualization.scenePath(imageFile, sceneDirectory);
  parameter Boolean saved(fixed = false);
  parameter Real bordered[pixels + 1, pixels + 1](each fixed = false) "Row 1: v in arcseconds, column 1: u in arcseconds, the rest the image";
initial algorithm
  // 1. An artificial observation of the true lens: rendered, blurred by the seeing, with Gaussian noise.
  noiseless := GeneralRelativity.Lensing.Inversion.renderImage(u, v, trueLens, trueSources, pointMass, 2, psfArcsec * arcsec);
  noiseSigma := noiseFraction * max(noiseless);
  observed := GeneralRelativity.Lensing.Inversion.syntheticImage(u, v, trueLens, trueSources, pointMass, noiseSigma, seed, 2, psfArcsec * arcsec);
  // 2. Where is the lens, how big is the ring, and what mass does that radius imply?
  // The bright pixels are fitted by a circle: its radius is the Einstein radius and its centre the lens (plus the source offset for arcs).
  (thetaE0, centre, ringSpread) := GeneralRelativity.Lensing.Inversion.einsteinRadiusFromImage(u, v, observed, threshold);
  massInitialSolar := GeneralRelativity.Lensing.massInEinsteinRadius(thetaE0, Dl, Ds, Dls) / GeneralRelativity.Constants.GM_sun;
  // 3. Simulate that estimate: a round lens of that radius at the centroid, a round exponential source behind it.
  lens0 := {thetaE0, centre[1], centre[2], 0.9, 0.0, 0.0, 0.0};
  sources0 := [centre[1], centre[2], max(observed), 0.3 * thetaE0, 0.9, 0.0, 1.0];
  firstGuess := GeneralRelativity.Lensing.Inversion.renderImage(u, v, lens0, sources0, pointMass, 2, psfArcsec * arcsec);
  // 4. How close is it?
  (sumSquaresInitial, chi2Initial, explainedInitial, correlationInitial) := GeneralRelativity.Lensing.Inversion.imageResidual(observed, firstGuess, noiseSigma);
  // 5. Iterate: Levenberg-Marquardt on every lens and source parameter until the simulated image matches.
  (lens, sources, chi2History, lensHistory, iterations, converged, lensSigma, sourcesSigma) := GeneralRelativity.Lensing.Inversion.fitLens(
    u, v, observed, noiseSigma, lens0, sources0, pointMass = pointMass, psfSigma = psfArcsec * arcsec, maxIterations = maxIterations);
  fitted := GeneralRelativity.Lensing.Inversion.renderImage(u, v, lens, sources, pointMass, 2, psfArcsec * arcsec);
  (sumSquaresFinal, chi2Final, explainedFinal, correlationFinal) := GeneralRelativity.Lensing.Inversion.imageResidual(observed, fitted, noiseSigma);
  massSolar := GeneralRelativity.Lensing.massInEinsteinRadius(lens[1], Dl, Ds, Dls) / GeneralRelativity.Constants.GM_sun;
  massSigmaSolar := 2.0 * lensSigma[1] / lens[1] * massSolar;
  massInApertureSolar := if pointMass then massSolar else massSolar * apertureM / (Dl * lens[1]);
  for k in 1:maxIterations + 1 loop
    massHistorySolar[k] := GeneralRelativity.Lensing.massInEinsteinRadius(lensHistory[k, 1], Dl, Ds, Dls) / GeneralRelativity.Constants.GM_sun;
  end for;
  massErrorSigmas := (massSolar - trueMassInRingSolar) / max(massSigmaSolar, 1.0e-300);
  // The scene: observed, first simulation, fitted simulation and residual images, the fit summary and the recovered lens.
  residual := (observed - fitted) / noiseSigma;
  initialDisplay := cat(1, lens0, sources0[1, :]) .* toDisplay;
  fittedDisplay := cat(1, lens, sources[1, :]) .* toDisplay;
  sigmaDisplay := cat(1, lensSigma, sourcesSigma[1, :]) .* toDisplay;
  truthDisplay := cat(1, trueLens, trueSources[1, :]) .* toDisplay;
  GeneralRelativity.Visualization.beginScene(path, "Lens inversion", "cartesian", 1.0, "arcsec", 1.0);
  GeneralRelativity.Visualization.writeView(path, "lensfit");
  GeneralRelativity.Visualization.writeLensImage(path, "observed", u / arcsec, v / arcsec, observed);
  GeneralRelativity.Visualization.writeLensImage(path, "initial", u / arcsec, v / arcsec, firstGuess);
  GeneralRelativity.Visualization.writeLensImage(path, "model", u / arcsec, v / arcsec, fitted);
  GeneralRelativity.Visualization.writeLensImage(path, "residual", u / arcsec, v / arcsec, residual, "sigma");
  GeneralRelativity.Visualization.writeLensFit(path, if pointMass then "point" else "galaxy",
    {"theta_E (arcsec)", "lens x (arcsec)", "lens y (arcsec)", "axis ratio q", "position angle (deg)", "shear", "shear angle (deg)",
     "source x (arcsec)", "source y (arcsec)", "source I_e", "source r_e (arcsec)", "source q", "source angle (deg)", "Sersic n"},
    initialDisplay, fittedDisplay, sigmaDisplay, truthDisplay, chi2History, massHistorySolar, massSolar, massSigmaSolar, massInitialSolar,
    trueMassInRingSolar, lens[1] / arcsec, {Dl, Ds, Dls} / GeneralRelativity.Constants.Mpc, iterations, converged, pixels^2, explainedFinal, correlationFinal);
  // The observation as one MATLAB v4 matrix (the v4 writer cannot append): u down the first column, v along the first row.
  bordered[1, 1] := 0.0;
  bordered[1, 2:pixels + 1] := v / arcsec;
  bordered[2:pixels + 1, 1] := u / arcsec;
  bordered[2:pixels + 1, 2:pixels + 1] := observed;
  saved := imageFile == "" or Modelica.Utilities.Streams.writeRealMatrix(imagePath, "lensImage", bordered, false, "4");
  GeneralRelativity.Visualization.writeLens(path, if pointMass then "point" else "galaxy", massInApertureSolar, apertureKpc, zLens, zSource, H0, OmegaM,
    lens[4], lens[5] / deg, lens[6], lens[7] / deg, lens[1] / arcsec, {Dl, Ds, Dls} / GeneralRelativity.Constants.Mpc);
equation
  assert(not checkRecovery or thetaE0 > 0.5 * trueThetaE and thetaE0 < 1.5 * trueThetaE, "The ring radius read off the image is not within 50% of the true Einstein radius");
  assert(not checkRecovery or chi2Final < chi2Initial, "Iterating did not improve on the first simulation");
  assert(not checkRecovery or chi2Final < 1.5 * pixels^2 and chi2Final > 0.7 * pixels^2, "The fitted simulation does not match the image to within the noise");
  assert(not checkRecovery or abs(lens[1] / trueThetaE - 1.0) < 0.02, "The fitted Einstein radius is not within 2% of the truth");
  assert(not checkRecovery or abs(massSolar - trueMassInRingSolar) < 3.0 * massSigmaSolar + 0.02 * trueMassInRingSolar, "The fitted mass is not within 3 standard errors (+2%) of the truth");
  assert(not checkRecovery or abs(massInApertureSolar / trueMassSolar - 1.0) < 0.04, "The mass within the aperture is not within 4% of the truth");
  annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 1.0e-8, Interval = 0.1),
    Documentation(info="<html>
<p>The inverse of <a href=\"modelica://GeneralRelativity.Examples.GalaxyLensing\">GalaxyLensing</a>. A galaxy of known mass
(the defaults are its 5&times;10<sup>11</sup> M<sub>&#9737;</sub> within 10 kpc at <code>z = 0.5 &rarr; 2</code>,
&theta;<sub>E</sub> = 1.25 arcsec) lenses an exponential source into arcs; the image is blurred by 0.06 arcsec seeing and 2%
noise is added. The model then forgets the truth and runs the steps of
<a href=\"modelica://GeneralRelativity.Lensing.Inversion\">Lensing.Inversion</a>:</p>
<ol>
<li>the image is the input (a file would do the same, see
<a href=\"modelica://GeneralRelativity.Examples.LensInversionFromFile\">LensInversionFromFile</a>);</li>
<li>the brightness centroid locates the lens and the mean radius of the bright pixels gives
<code>&theta;<sub>E,0</sub></code>, hence the projected mass inside it, <code>&pi;(D<sub>l</sub>&theta;<sub>E</sub>)&sup2;&Sigma;<sub>cr</sub></code>
(<code>massInitialSolar</code>);</li>
<li>a round lens of that radius with a round source behind it is simulated (<code>firstGuess</code>);</li>
<li>&chi;&sup2;, explained fraction and correlation against the observation (<code>chi2Initial</code>, ...);</li>
<li>Levenberg-Marquardt refines all 14 parameters until the simulation matches (<code>lens</code>, <code>sources</code>,
<code>chi2History</code>); <code>massSolar &plusmn; massSigmaSolar</code> is the answer and <code>massErrorSigmas</code>
grades it against the truth.</li>
</ol>
<p>The assertions require the ring estimate within 50%, the fitted Einstein radius within 2%, the mass within three standard
errors and the final &chi;&sup2; consistent with the noise. The scene opens in the viewer's <b>lensfit</b> style (observed,
first and fitted simulations, residual in units of &sigma;, the &chi;&sup2; and mass history against the truth) and carries a
<code>lens</code> record, so the <b>lensing</b> style shows a synthetic deep field through the recovered lens.</p>
</html>"));
end LensInversion;
