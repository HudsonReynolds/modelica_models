within GeneralRelativity.Examples;
model LensInversionFromFile "Lens inversion of an image read from a file: the mass that must be distorting the light, with its standard error, and the simulated image that matches"
  parameter String imageFile = GeneralRelativity.Visualization.scenePath("lens_inversion_observed.mat")
    "MATLAB v4 file with the bordered matrix lensImage (first column u and first row v in arcseconds, the rest the image): written by LensInversion, or by tools/lens_image_to_mat.py from a picture";
  parameter Real zLens = 0.5 "Lens redshift";
  parameter Real zSource = 2.0 "Source redshift";
  parameter Real H0 = 70.0 "Hubble constant, km/s/Mpc";
  parameter Real OmegaM = 0.3 "Matter density (flat universe)";
  parameter Boolean pointMass = false "true: point-mass lens; false: isothermal ellipse with external shear";
  parameter Real apertureKpc = 10.0 "Aperture radius for reporting a galaxy mass, as in GalaxyLensing";
  parameter Real psfArcsec = 0.06 "Gaussian point-spread width of the picture";
  parameter Real noiseSigmaInput = 0.0 "Pixel noise; 0 estimates it from the pixels below the threshold";
  parameter Real threshold = 0.2 "Fraction of the peak below which pixels are background";
  parameter Integer maxIterations = 30;
  parameter String sceneFile = "lens_inversion_file.ndjson" "Scene file name";
  parameter String sceneDirectory = "" "Output folder; empty uses GR_SCENE_DIR, else ~/Documents/Visualization";
  final parameter Integer borderedSize[2] = Modelica.Utilities.Streams.readMatrixSize(imageFile, "lensImage");
  final parameter Integer dims[2] = borderedSize - {1, 1} "Image size";
  final parameter Real bordered[dims[1] + 1, dims[2] + 1] = Modelica.Utilities.Streams.readRealMatrix(imageFile, "lensImage", dims[1] + 1, dims[2] + 1);
  final parameter Real observed[dims[1], dims[2]] = bordered[2:dims[1] + 1, 2:dims[2] + 1] "The input image";
  final parameter Real arcsec = GeneralRelativity.Constants.arcsec;
  final parameter Real deg = GeneralRelativity.Constants.pi / 180.0;
  final parameter Real u[dims[1]](each unit = "rad") = bordered[2:dims[1] + 1, 1] * arcsec;
  final parameter Real v[dims[2]](each unit = "rad") = bordered[1, 2:dims[2] + 1] * arcsec;
  final parameter Real Dl(unit = "m") = GeneralRelativity.Lensing.angularDiameterDistance(0.0, zLens, H0, OmegaM);
  final parameter Real Ds(unit = "m") = GeneralRelativity.Lensing.angularDiameterDistance(0.0, zSource, H0, OmegaM);
  final parameter Real Dls(unit = "m") = GeneralRelativity.Lensing.angularDiameterDistance(zLens, zSource, H0, OmegaM);
  final parameter Real apertureM(unit = "m") = apertureKpc * 1.0e-3 * GeneralRelativity.Constants.Mpc;
  final parameter Real toDisplay[14] = {1.0 / arcsec, 1.0 / arcsec, 1.0 / arcsec, 1.0, 1.0 / deg, 1.0, 1.0 / deg, 1.0 / arcsec, 1.0 / arcsec, 1.0, 1.0 / arcsec, 1.0, 1.0 / deg, 1.0};
  parameter Real noiseSigma(fixed = false) "Pixel noise used for chi2";
  parameter Real centre[2](each fixed = false, each unit = "rad");
  parameter Real thetaE0(fixed = false, unit = "rad") "Einstein radius read off the image";
  parameter Real ringSpread(fixed = false, unit = "rad");
  parameter Real massInitialSolar(fixed = false) "Radius-based mass estimate inside thetaE0, solar masses";
  parameter Real lens0[7](each fixed = false);
  parameter Real sources0[1, 7](each fixed = false);
  parameter Real firstGuess[dims[1], dims[2]](each fixed = false);
  parameter Real chi2Initial(fixed = false);
  parameter Real lens[7](each fixed = false) "Fitted lens";
  parameter Real sources[1, 7](each fixed = false);
  parameter Real chi2History[maxIterations + 1](each fixed = false);
  parameter Real lensHistory[maxIterations + 1, 7](each fixed = false);
  parameter Integer iterations(fixed = false);
  parameter Boolean converged(fixed = false);
  parameter Real lensSigma[7](each fixed = false);
  parameter Real sourcesSigma[1, 7](each fixed = false);
  parameter Real fitted[dims[1], dims[2]](each fixed = false);
  parameter Real chi2Final(fixed = false);
  parameter Real explainedFinal(fixed = false);
  parameter Real correlationFinal(fixed = false);
  parameter Real massSolar(fixed = false) "Projected mass inside the Einstein radius, solar masses";
  parameter Real massSigmaSolar(fixed = false);
  parameter Real massInApertureSolar(fixed = false);
  parameter Real massHistorySolar[maxIterations + 1](each fixed = false);
protected
  parameter Real sumSquares(fixed = false);
  parameter Real explainedInitial(fixed = false);
  parameter Real correlationInitial(fixed = false);
  parameter Real residual[dims[1], dims[2]](each fixed = false);
  parameter Real initialDisplay[14](each fixed = false);
  parameter Real fittedDisplay[14](each fixed = false);
  parameter Real sigmaDisplay[14](each fixed = false);
  parameter Real peak(fixed = false);
  parameter Real backgroundSum(fixed = false);
  parameter Real backgroundSum2(fixed = false);
  parameter Integer backgroundCount(fixed = false);
  parameter String path = GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory);
initial algorithm
  // 1. Noise level: the scatter of the pixels that hold no arc light.
  peak := max(observed);
  backgroundSum := 0.0;
  backgroundSum2 := 0.0;
  backgroundCount := 0;
  for i in 1:dims[1] loop
    for j in 1:dims[2] loop
      if observed[i, j] < threshold * peak then
        backgroundSum := backgroundSum + observed[i, j];
        backgroundSum2 := backgroundSum2 + observed[i, j]^2;
        backgroundCount := backgroundCount + 1;
      end if;
    end for;
  end for;
  noiseSigma := if noiseSigmaInput > 0.0 then noiseSigmaInput
    else max(sqrt(max(backgroundSum2 / max(backgroundCount, 1) - (backgroundSum / max(backgroundCount, 1))^2, 0.0)), 1.0e-3 * peak);
  // 2. Ring radius and the mass it implies.
  // The bright pixels are fitted by a circle: its radius is the Einstein radius and its centre the lens (plus the source offset for arcs).
  (thetaE0, centre, ringSpread) := GeneralRelativity.Lensing.Inversion.einsteinRadiusFromImage(u, v, observed, threshold);
  massInitialSolar := GeneralRelativity.Lensing.massInEinsteinRadius(thetaE0, Dl, Ds, Dls) / GeneralRelativity.Constants.GM_sun;
  // 3 and 4. Simulate the estimate and compare.
  lens0 := {thetaE0, centre[1], centre[2], 0.9, 0.0, 0.0, 0.0};
  sources0 := [centre[1], centre[2], peak, 0.3 * thetaE0, 0.9, 0.0, 1.0];
  firstGuess := GeneralRelativity.Lensing.Inversion.renderImage(u, v, lens0, sources0, pointMass, 2, psfArcsec * arcsec);
  (sumSquares, chi2Initial, explainedInitial, correlationInitial) := GeneralRelativity.Lensing.Inversion.imageResidual(observed, firstGuess, noiseSigma);
  // 5. Iterate.
  (lens, sources, chi2History, lensHistory, iterations, converged, lensSigma, sourcesSigma) := GeneralRelativity.Lensing.Inversion.fitLens(
    u, v, observed, noiseSigma, lens0, sources0, pointMass = pointMass, psfSigma = psfArcsec * arcsec, maxIterations = maxIterations);
  fitted := GeneralRelativity.Lensing.Inversion.renderImage(u, v, lens, sources, pointMass, 2, psfArcsec * arcsec);
  (sumSquares, chi2Final, explainedFinal, correlationFinal) := GeneralRelativity.Lensing.Inversion.imageResidual(observed, fitted, noiseSigma);
  massSolar := GeneralRelativity.Lensing.massInEinsteinRadius(lens[1], Dl, Ds, Dls) / GeneralRelativity.Constants.GM_sun;
  massSigmaSolar := 2.0 * lensSigma[1] / lens[1] * massSolar;
  massInApertureSolar := if pointMass then massSolar else massSolar * apertureM / (Dl * lens[1]);
  for k in 1:maxIterations + 1 loop
    massHistorySolar[k] := GeneralRelativity.Lensing.massInEinsteinRadius(lensHistory[k, 1], Dl, Ds, Dls) / GeneralRelativity.Constants.GM_sun;
  end for;
  residual := (observed - fitted) / noiseSigma;
  initialDisplay := cat(1, lens0, sources0[1, :]) .* toDisplay;
  fittedDisplay := cat(1, lens, sources[1, :]) .* toDisplay;
  sigmaDisplay := cat(1, lensSigma, sourcesSigma[1, :]) .* toDisplay;
  GeneralRelativity.Visualization.beginScene(path, "Lens inversion of " + imageFile, "cartesian", 1.0, "arcsec", 1.0);
  GeneralRelativity.Visualization.writeView(path, "lensfit");
  GeneralRelativity.Visualization.writeLensImage(path, "observed", u / arcsec, v / arcsec, observed);
  GeneralRelativity.Visualization.writeLensImage(path, "initial", u / arcsec, v / arcsec, firstGuess);
  GeneralRelativity.Visualization.writeLensImage(path, "model", u / arcsec, v / arcsec, fitted);
  GeneralRelativity.Visualization.writeLensImage(path, "residual", u / arcsec, v / arcsec, residual, "sigma");
  GeneralRelativity.Visualization.writeLensFit(path, if pointMass then "point" else "galaxy",
    {"theta_E (arcsec)", "lens x (arcsec)", "lens y (arcsec)", "axis ratio q", "position angle (deg)", "shear", "shear angle (deg)",
     "source x (arcsec)", "source y (arcsec)", "source I_e", "source r_e (arcsec)", "source q", "source angle (deg)", "Sersic n"},
    initialDisplay, fittedDisplay, sigmaDisplay, fill(0.0, 0), chi2History, massHistorySolar, massSolar, massSigmaSolar, massInitialSolar,
    0.0, lens[1] / arcsec, {Dl, Ds, Dls} / GeneralRelativity.Constants.Mpc, iterations, converged, dims[1] * dims[2], explainedFinal, correlationFinal);
  GeneralRelativity.Visualization.writeLens(path, if pointMass then "point" else "galaxy", massInApertureSolar, apertureKpc, zLens, zSource, H0, OmegaM,
    lens[4], lens[5] / deg, lens[6], lens[7] / deg, lens[1] / arcsec, {Dl, Ds, Dls} / GeneralRelativity.Constants.Mpc);
equation
  assert(chi2Final <= chi2Initial, "Iterating did not improve on the first simulation");
  assert(lens[1] > 0.0 and massSolar > 0.0, "The fit found no lens");
  annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 1.0e-8, Interval = 0.1),
    Documentation(info="<html>
<p>The same five steps as <a href=\"modelica://GeneralRelativity.Examples.LensInversion\">LensInversion</a> for an image whose
true lens is unknown. The file holds one MATLAB v4 matrix, <code>lensImage</code>: its first column (from row 2) is <code>u</code> and
its first row (from column 2) is <code>v</code>, the pixel-centre angles in arcseconds, and the rest is the surface brightness with the
lens light subtracted. <code>LensInversion</code> writes such a file of
its synthetic observation, so by default this model re-analyses it; <code>tools/lens_image_to_mat.py</code> makes one from a PNG or
JPEG given the pixel scale. The redshifts and seeing must be supplied: they are not in a picture. The noise is estimated from the
pixels below the threshold unless given.</p>
<p>Results: <code>massInitialSolar</code> from the ring radius, then <code>massSolar &plusmn; massSigmaSolar</code> (projected mass
inside the Einstein radius) and <code>massInApertureSolar</code> (an isothermal galaxy's mass within the aperture) from the fit,
with the parameter standard errors in <code>lensSigma</code> and <code>sourcesSigma</code>, <code>chi2Final</code> against
the number of pixels, and the scene for the viewer's <b>lensfit</b> style.</p>
</html>"));
end LensInversionFromFile;
