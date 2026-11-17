within GeneralRelativity.Examples;
model LensInversion "Lens inversion on a synthetic image: read the Einstein radius off the arcs, estimate the mass, simulate, compare and iterate until the simulation matches; the true mass is known and the recovery is asserted"
  // The truth, used only to make the image and to grade the answer.
  parameter Real trueMassSolar = 5.0e11 "Mass that sets the deflection scale b (as in GalaxyLensing: of the isothermal sphere with b, within the aperture; total for a point)";
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
  parameter Integer maxIterations = 30;
  parameter Boolean checkRecovery = true "false switches off the assertions that grade the fit against the truth";
  parameter String sceneFile = "lens_inversion.ndjson" "Scene file name";
  parameter String sceneDirectory = "" "Output folder; empty uses GR_SCENE_DIR, else ~/Documents/Visualization";
  parameter String imageFile = "lens_inversion_observed.mat" "The observed image is also saved here as the bordered matrix lensImage (MATLAB v6), the input of LensInversionFromFile; empty to skip";
  final parameter Real arcsec = GeneralRelativity.Constants.arcsec;
  final parameter Real deg = GeneralRelativity.Constants.pi / 180.0;
  final parameter Real Dl(unit = "m") = GeneralRelativity.Lensing.angularDiameterDistance(0.0, zLens, H0, OmegaM);
  final parameter Real Ds(unit = "m") = GeneralRelativity.Lensing.angularDiameterDistance(0.0, zSource, H0, OmegaM);
  final parameter Real Dls(unit = "m") = GeneralRelativity.Lensing.angularDiameterDistance(zLens, zSource, H0, OmegaM);
  final parameter Real apertureM(unit = "m") = apertureKpc * 1.0e-3 * GeneralRelativity.Constants.Mpc;
  final parameter Real trueGM(unit = "m3/s2") = trueMassSolar * GeneralRelativity.Constants.GM_sun;
  final parameter Real trueB(unit = "rad") = if pointMass then GeneralRelativity.Lensing.einsteinRadiusPointMass(trueGM, Dl, Ds, Dls)
    else GeneralRelativity.Lensing.einsteinRadiusSIS(sqrt(trueGM / (GeneralRelativity.Constants.pi * apertureM)), Ds, Dls) "True deflection scale";
  final parameter Real trueLens[7] = {trueB, lensCentreArcsec[1] * arcsec, lensCentreArcsec[2] * arcsec, axisRatio, positionAngleDeg * deg, shear, shearAngleDeg * deg};
  final parameter Real trueSource[7] = {sourceArcsec[1] * arcsec, sourceArcsec[2] * arcsec, 1.0, sourceRadiusArcsec * arcsec, sourceAxisRatio, sourceAngleDeg * deg, sersicIndex};
  final parameter Real u[pixels](each unit = "rad") = {(-fieldArcsec + (i - 0.5) * 2.0 * fieldArcsec / pixels) * arcsec for i in 1:pixels} "Pixel centres";
  final parameter Real v[pixels](each unit = "rad") = u;
  // Step 1: the image.
  parameter Real noiseSigma(fixed = false) "Pixel noise standard deviation";
  parameter Real observed[pixels, pixels](each fixed = false) "The input image";
  // Steps 2 to 5.
  parameter Real ringRadiusArcsec(fixed = false) "Einstein radius read off the image";
  parameter Real massInitialSolar(fixed = false) "Mass inside that radius: the radius-based estimate";
  parameter Real lens[7](each fixed = false) "Fitted lens";
  parameter Real sources[1, 7](each fixed = false) "Fitted source";
  parameter Real lensSigma[7](each fixed = false) "Standard errors of the lens parameters";
  parameter Real thetaEArcsec(fixed = false) "Fitted equal-area Einstein radius";
  parameter Real thetaESigmaArcsec(fixed = false);
  parameter Real massSolar(fixed = false) "Fitted projected mass inside the critical curve, solar masses";
  parameter Real massSigmaSolar(fixed = false) "Its standard error";
  parameter Real massInApertureSolar(fixed = false) "Mass within the aperture implied by the fitted profile";
  parameter Real trueThetaEArcsec(fixed = false) "Equal-area Einstein radius of the true lens";
  parameter Real trueMassInRingSolar(fixed = false) "Mass inside the true critical curve";
  parameter Real trueMassInApertureSolar(fixed = false) "True mass within the aperture";
  parameter Real chi2Initial(fixed = false) "Cost of the first simulation";
  parameter Real chi2Final(fixed = false) "Cost of the fit; about pixels^2 when the model and the noise are right";
  parameter Real explainedFinal(fixed = false);
  parameter Real correlationFinal(fixed = false);
  parameter Integer iterations(fixed = false);
  parameter Boolean converged(fixed = false);
  parameter Real massErrorSigmas(fixed = false) "(fitted - true mass inside the critical curve) / standard error";
protected
  parameter String path = GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory);
  parameter String imagePath = if imageFile == "" then "" else GeneralRelativity.Visualization.scenePath(imageFile, sceneDirectory);
  parameter Boolean saved(fixed = false);
  parameter Real bordered[pixels + 1, pixels + 1](each fixed = false) "Row 1: v in arcseconds, column 1: u in arcseconds, the rest the image";
  parameter Real noiseless[pixels, pixels](each fixed = false);
  parameter Real ringRadius(fixed = false);
  parameter Real thetaE(fixed = false);
  parameter Real thetaESigma(fixed = false);
  parameter Real trueThetaE(fixed = false);
initial algorithm
  // 1. An artificial observation of the true lens: rendered, blurred by the seeing, with Gaussian noise.
  noiseless := GeneralRelativity.Lensing.Inversion.renderImage(u, v, trueLens, {trueSource}, pointMass, 2, psfArcsec * arcsec);
  noiseSigma := noiseFraction * max(noiseless);
  observed := GeneralRelativity.Lensing.Inversion.syntheticImage(u, v, trueLens, {trueSource}, pointMass, noiseSigma, seed, 2, psfArcsec * arcsec);
  // 2 to 5, written to the scene.
  GeneralRelativity.Visualization.beginScene(path, "Lens inversion", "cartesian", 1.0, "arcsec", 1.0);
  GeneralRelativity.Visualization.writeView(path, "lensfit");
  (ringRadius, massInitialSolar, lens, sources, lensSigma, thetaE, thetaESigma, massSolar, massSigmaSolar, massInApertureSolar, trueThetaE, trueMassInRingSolar,
    chi2Initial, chi2Final, explainedFinal, correlationFinal, iterations, converged) := GeneralRelativity.Visualization.writeLensInversion(path, "", "Synthetic galaxy lens",
    "An isothermal ellipse with shear lenses an exponential source; seeing and 2% noise are added, then the lens is recovered from the pixels alone.",
    u, v, observed, noiseSigma, pointMass, psfArcsec * arcsec, {Dl, Ds, Dls}, apertureM, trueLens, trueSource, maxIterations, 4, 0.2, true, zLens, zSource, H0, OmegaM);
  ringRadiusArcsec := ringRadius / arcsec;
  thetaEArcsec := thetaE / arcsec;
  thetaESigmaArcsec := thetaESigma / arcsec;
  trueThetaEArcsec := trueThetaE / arcsec;
  trueMassInApertureSolar := if pointMass then trueMassInRingSolar else trueMassInRingSolar * apertureM / (Dl * trueThetaE);
  massErrorSigmas := (massSolar - trueMassInRingSolar) / max(massSigmaSolar, 1.0e-300);
  // The observation as one MATLAB matrix: u down the first column, v along the first row. Format v6: the v4 writer of
  // OpenModelica 1.28 on macOS writes an empty file (readRealMatrix reads v4 and v6 alike).
  bordered[1, 1] := 0.0;
  bordered[1, 2:pixels + 1] := v / arcsec;
  bordered[2:pixels + 1, 1] := u / arcsec;
  bordered[2:pixels + 1, 2:pixels + 1] := observed;
  saved := imageFile == "" or Modelica.Utilities.Streams.writeRealMatrix(imagePath, "lensImage", bordered, false, "6");
equation
  assert(not checkRecovery or abs(ringRadius / trueThetaE - 1.0) < 0.1, "The ring radius read off the image is not within 10% of the true Einstein radius");
  assert(not checkRecovery or chi2Final < chi2Initial, "Iterating did not improve on the first simulation");
  assert(not checkRecovery or (chi2Final < 1.5 * pixels^2 and chi2Final > 0.7 * pixels^2), "The fitted simulation does not match the image to within the noise");
  assert(not checkRecovery or abs(lens[1] / trueB - 1.0) < 0.02, "The fitted deflection scale is not within 2% of the truth");
  assert(not checkRecovery or abs(massSolar - trueMassInRingSolar) < 3.0 * massSigmaSolar + 0.02 * trueMassInRingSolar, "The fitted mass is not within 3 standard errors (+2%) of the truth");
  assert(not checkRecovery or abs(massInApertureSolar / trueMassInApertureSolar - 1.0) < 0.04, "The mass within the aperture is not within 4% of the truth");
  assert(saved, "The observed image could not be saved");
  annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 1.0e-8, Interval = 0.1),
    Documentation(info="<html>
<p>The inverse of <a href=\"modelica://GeneralRelativity.Examples.GalaxyLensing\">GalaxyLensing</a>. A galaxy with GalaxyLensing's
defaults (deflection scale of 5&times;10<sup>11</sup> M<sub>&#9737;</sub> within 10 kpc at <code>z = 0.5 &rarr; 2</code>, axis ratio 0.75)
lenses an exponential source into arcs; the image is blurred by 0.06 arcsec seeing and 2% noise is added. The model then forgets
the truth and runs the steps of <a href=\"modelica://GeneralRelativity.Lensing.Inversion\">Lensing.Inversion</a> through
<a href=\"modelica://GeneralRelativity.Visualization.writeLensInversion\">Visualization.writeLensInversion</a>:</p>
<ol>
<li>the image is the input (also saved for <a href=\"modelica://GeneralRelativity.Examples.LensInversionFromFile\">LensInversionFromFile</a>);</li>
<li>the circle through the arcs gives the Einstein radius <code>ringRadiusArcsec</code>, hence the mass inside it
(<code>massInitialSolar</code>);</li>
<li>a round lens of that radius is simulated;</li>
<li>&chi;&sup2; against the observation (<code>chi2Initial</code>);</li>
<li>short fits from two candidate lens positions, then Levenberg-Marquardt on all 14 parameters until the simulation matches:
<code>massSolar &plusmn; massSigmaSolar</code> inside the critical curve, graded by <code>massErrorSigmas</code>.</li>
</ol>
<p>The mass inside the critical curve of an isothermal ellipse is <code>&pi;qb&sup2;&Sigma;<sub>cr</sub>D<sub>l</sub>&sup2;</code>
(<a href=\"modelica://GeneralRelativity.Lensing.Inversion.einsteinMass\">einsteinMass</a>), so with q = 0.75 it is three quarters of
what the round sphere with the same <code>b</code> would hold. The scene opens in the viewer's <b>lensfit</b> style and carries a
<code>lens</code> record for the <b>lensing</b> style.</p>
</html>"));
end LensInversion;
