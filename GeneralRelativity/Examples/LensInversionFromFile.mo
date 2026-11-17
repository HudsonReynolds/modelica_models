within GeneralRelativity.Examples;
model LensInversionFromFile "Lens inversion of an image read from a file: the mass that must be distorting the light, with its standard error, and the simulated image that matches"
  parameter String imageFile = GeneralRelativity.Visualization.scenePath("lens_inversion_observed.mat")
    "MATLAB file (v4 or v6) with the bordered matrix lensImage (first column u and first row v in arcseconds, the rest the image): written by LensInversion, or by tools/lens_image_to_mat.py from a picture";
  parameter String title = "Image from a file" "Title shown by the viewer";
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
  // A missing file would crash the compiler (the size is read at translation): fall back to a dummy 2 x 2 image
  // and stop at initialization with a message instead.
  final parameter Boolean imageFound = Modelica.Utilities.Files.exist(imageFile) "The image file exists";
  final parameter Integer borderedSize[2] = if imageFound then Modelica.Utilities.Streams.readMatrixSize(imageFile, "lensImage") else {3, 3};
  final parameter Integer dims[2] = borderedSize - {1, 1} "Image size";
  final parameter Real bordered[dims[1] + 1, dims[2] + 1] = if imageFound then Modelica.Utilities.Streams.readRealMatrix(imageFile, "lensImage", dims[1] + 1, dims[2] + 1)
    else [0, 1, 2; 1, 0, 0; 2, 0, 0];
  final parameter Real observed[dims[1], dims[2]] = bordered[2:dims[1] + 1, 2:dims[2] + 1] "The input image";
  final parameter Real arcsec = GeneralRelativity.Constants.arcsec;
  final parameter Real u[dims[1]](each unit = "rad") = bordered[2:dims[1] + 1, 1] * arcsec;
  final parameter Real v[dims[2]](each unit = "rad") = bordered[1, 2:dims[2] + 1] * arcsec;
  final parameter Real Dl(unit = "m") = GeneralRelativity.Lensing.angularDiameterDistance(0.0, zLens, H0, OmegaM);
  final parameter Real Ds(unit = "m") = GeneralRelativity.Lensing.angularDiameterDistance(0.0, zSource, H0, OmegaM);
  final parameter Real Dls(unit = "m") = GeneralRelativity.Lensing.angularDiameterDistance(zLens, zSource, H0, OmegaM);
  final parameter Real apertureM(unit = "m") = apertureKpc * 1.0e-3 * GeneralRelativity.Constants.Mpc;
  parameter Real noiseSigma(fixed = false) "Pixel noise used for chi2";
  parameter Real ringRadiusArcsec(fixed = false) "Einstein radius read off the image";
  parameter Real massInitialSolar(fixed = false) "Radius-based mass estimate, solar masses";
  parameter Real lens[7](each fixed = false) "Fitted lens";
  parameter Real sources[1, 7](each fixed = false);
  parameter Real lensSigma[7](each fixed = false);
  parameter Real thetaEArcsec(fixed = false) "Fitted equal-area Einstein radius";
  parameter Real thetaESigmaArcsec(fixed = false);
  parameter Real massSolar(fixed = false) "Projected mass inside the critical curve, solar masses";
  parameter Real massSigmaSolar(fixed = false);
  parameter Real massInApertureSolar(fixed = false);
  parameter Real chi2Initial(fixed = false);
  parameter Real chi2Final(fixed = false) "About the number of pixels when the fit reaches the noise";
  parameter Real explainedFinal(fixed = false);
  parameter Real correlationFinal(fixed = false);
  parameter Integer iterations(fixed = false);
  parameter Boolean converged(fixed = false);
protected
  parameter Real ringRadius(fixed = false);
  parameter Real thetaE(fixed = false);
  parameter Real thetaESigma(fixed = false);
  parameter Real unknownThetaE(fixed = false);
  parameter Real unknownMass(fixed = false);
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
  // 2 to 5, written to the scene.
  GeneralRelativity.Visualization.beginScene(path, "Lens inversion", "cartesian", 1.0, "arcsec", 1.0);
  GeneralRelativity.Visualization.writeView(path, "lensfit");
  (ringRadius, massInitialSolar, lens, sources, lensSigma, thetaE, thetaESigma, massSolar, massSigmaSolar, massInApertureSolar, unknownThetaE, unknownMass,
    chi2Initial, chi2Final, explainedFinal, correlationFinal, iterations, converged) := GeneralRelativity.Visualization.writeLensInversion(path, "", title,
    "Read from " + imageFile + "; the true lens is unknown.", u, v, observed, noiseSigma, pointMass, psfArcsec * arcsec, {Dl, Ds, Dls}, apertureM,
    fill(0.0, 0), fill(0.0, 0), maxIterations, 4, threshold, true, zLens, zSource, H0, OmegaM);
  ringRadiusArcsec := ringRadius / arcsec;
  thetaEArcsec := thetaE / arcsec;
  thetaESigmaArcsec := thetaESigma / arcsec;
initial equation
  assert(imageFound, "No image file " + imageFile + ": run Examples.LensInversion first (it writes one), or make one from a picture with tools/lens_image_to_mat.py");
equation
  assert(chi2Final <= chi2Initial, "Iterating did not improve on the first simulation");
  assert(lens[1] > 0.0 and massSolar > 0.0, "The fit found no lens");
  annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 1.0e-8, Interval = 0.1),
    Documentation(info="<html>
<p>The same steps as <a href=\"modelica://GeneralRelativity.Examples.LensInversion\">LensInversion</a> for an image whose
true lens is unknown. The file holds one MATLAB matrix (v4 or v6), <code>lensImage</code>: its first column (from row 2) is <code>u</code> and
its first row (from column 2) is <code>v</code>, the pixel-centre angles in arcseconds, and the rest is the surface brightness with the
lens light subtracted. <code>LensInversion</code> writes such a file of its synthetic observation, so by default this model
re-analyses it; <code>tools/lens_image_to_mat.py</code> makes one from a PNG or JPEG given the pixel scale. The redshifts and
seeing must be supplied: they are not in a picture. The noise is estimated from the pixels below the threshold unless given.</p>
<p>Results: <code>massInitialSolar</code> from the ring radius, then <code>massSolar &plusmn; massSigmaSolar</code> (projected mass
inside the critical curve) and <code>massInApertureSolar</code> (an isothermal galaxy's mass within the aperture) from the fit,
<code>chi2Final</code> against the number of pixels, and the scene for the viewer's <b>lensfit</b> style.</p>
</html>"));
end LensInversionFromFile;
