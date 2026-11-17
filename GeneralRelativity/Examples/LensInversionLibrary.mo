within GeneralRelativity.Examples;
model LensInversionLibrary "A library of synthetic lensed images (ring, double, giant arc, Einstein cross, point mass, noisy, group-scale, low-mass), each inverted to the mass behind it and written to one scene with a case gallery"
  constant Integer nCases = 8;
  parameter String caseNames[nCases] = {"ring", "double", "arc", "cross", "point", "noisy", "group", "small"} "Short identifiers, used as the viewer's case keys";
  parameter String titles[nCases] = {
    "Einstein ring",
    "Double image",
    "Giant arc",
    "Einstein cross",
    "Point-mass lens",
    "Faint, blurred ring",
    "Group-scale lens",
    "Low-mass galaxy"};
  parameter String descriptions[nCases] = {
    "A source almost directly behind a round galaxy is stretched into a nearly complete ring at the Einstein radius.",
    "A source well off the axis gives a bright outer image and a faint inner counter-image on the opposite side of the lens.",
    "A flattened galaxy with external shear and a source near a cusp of its caustic: three images merge into one long tangential arc.",
    "A compact source, like a quasar, close to the centre of an elongated lens: four images in a cross around the galaxy.",
    "A compact mass such as a black hole or a star cluster: two images on opposite sides of the lens, the inner one fainter and mirror-reversed.",
    "A lighter galaxy seen through twice the seeing and with four times the noise: the mass is still recovered, with a larger error bar.",
    "A group of galaxies at low redshift lensing a large galaxy at z = 3: arcs three arcseconds from the centre.",
    "A small galaxy with an Einstein radius of a quarter of an arcsecond, observed at high resolution."};
  parameter Boolean pointMass[nCases] = {false, false, false, false, true, false, false, false} "true: point-mass lens";
  parameter Real massSolar[nCases] = {5.0e11, 4.0e11, 6.0e11, 5.0e11, 2.0e11, 3.0e11, 1.0e12, 1.2e11}
    "Mass that sets the deflection scale b, as in GalaxyLensing: of the isothermal sphere with b within the aperture, or the total of a point";
  parameter Real zLens[nCases] = {0.3, 0.5, 0.6, 0.4, 0.5, 0.7, 0.35, 0.4};
  parameter Real zSource[nCases] = {1.5, 2.0, 2.5, 1.8, 2.0, 2.2, 3.0, 1.5};
  parameter Real axisRatio[nCases] = {0.9, 0.85, 0.55, 0.6, 1.0, 0.8, 0.8, 0.75} "Lens axis ratio (ignored by a point mass)";
  parameter Real positionAngleDeg[nCases] = {20.0, 70.0, 35.0, 45.0, 0.0, 110.0, 150.0, 60.0};
  parameter Real shear[nCases] = {0.02, 0.03, 0.05, 0.06, 0.0, 0.03, 0.1, 0.02};
  parameter Real shearAngleDeg[nCases] = {60.0, 10.0, 100.0, 45.0, 0.0, 30.0, 70.0, 120.0};
  parameter Real lensCentreArcsec[nCases, 2] = [0.0, 0.0; -0.1, 0.05; 0.05, -0.05; 0.0, 0.0; 0.1, 0.0; -0.05, 0.08; 0.2, -0.1; 0.0, 0.02];
  parameter Real sourceArcsec[nCases, 2] = [0.04, -0.03; 0.6, 0.3; 0.21, 0.13; 0.04, 0.02; 0.35, -0.2; 0.15, 0.1; 0.5, 0.4; 0.05, 0.08]
    "Source position in the image frame";
  parameter Real sourceRadiusArcsec[nCases] = {0.2, 0.12, 0.1, 0.02, 0.05, 0.15, 0.3, 0.06};
  parameter Real sourceAxisRatio[nCases] = {0.8, 0.7, 0.6, 0.9, 0.9, 0.7, 0.6, 0.8};
  parameter Real sourceAngleDeg[nCases] = {30.0, 120.0, 80.0, 10.0, 45.0, 160.0, 20.0, 100.0};
  parameter Real sersicIndex[nCases] = {1.0, 1.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.5};
  parameter Real fieldArcsec[nCases] = {3.0, 3.0, 3.5, 2.5, 2.5, 2.5, 5.0, 1.4} "Half-width of each image";
  parameter Real psfArcsec[nCases] = {0.06, 0.06, 0.06, 0.06, 0.05, 0.12, 0.1, 0.03} "Seeing (Gaussian sigma), assumed known";
  parameter Real noiseFraction[nCases] = {0.02, 0.02, 0.02, 0.02, 0.02, 0.08, 0.02, 0.02} "Pixel noise as a fraction of the brightest noiseless pixel";
  parameter Integer pixels = 64 "Image size per axis";
  parameter Real apertureKpc = 10.0 "Aperture radius that defines a galaxy mass";
  parameter Real H0 = 70.0 "Hubble constant, km/s/Mpc";
  parameter Real OmegaM = 0.3 "Matter density (flat universe)";
  parameter Integer maxIterations = 40 "Iterations of the final refinement of each case";
  parameter Boolean checkRecovery = true "false switches off the assertions that grade each fit against the truth";
  parameter String sceneFile = "lens_library.ndjson" "Scene file name";
  parameter String sceneDirectory = "" "Output folder; empty uses GR_SCENE_DIR, else ~/Documents/Visualization";
  // Results, one per case.
  parameter Real trueThetaEArcsec[nCases](each fixed = false) "Equal-area Einstein radius of the true lens";
  parameter Real trueMassInRingSolar[nCases](each fixed = false) "Mass inside the true critical curve";
  parameter Real ringRadiusArcsec[nCases](each fixed = false) "Einstein radius read off the image";
  parameter Real massInitialSolar[nCases](each fixed = false) "Mass implied by that radius";
  parameter Real thetaEArcsec[nCases](each fixed = false) "Fitted equal-area Einstein radius";
  parameter Real thetaESigmaArcsec[nCases](each fixed = false) "Its standard error";
  parameter Real massFitSolar[nCases](each fixed = false) "Fitted projected mass inside the critical curve";
  parameter Real massSigmaSolar[nCases](each fixed = false) "Its standard error";
  parameter Real massErrorSigmas[nCases](each fixed = false) "(fitted - true mass) / standard error";
  parameter Real chi2Initial[nCases](each fixed = false) "Cost of the first simulation";
  parameter Real chi2PerPixel[nCases](each fixed = false) "Final cost per pixel, about 1 when the fit reaches the noise";
  parameter Integer iterations[nCases](each fixed = false);
  parameter Boolean converged[nCases](each fixed = false);
protected
  final parameter Real arcsec = GeneralRelativity.Constants.arcsec;
  final parameter Real deg = GeneralRelativity.Constants.pi / 180.0;
  final parameter Real apertureM = apertureKpc * 1.0e-3 * GeneralRelativity.Constants.Mpc;
  parameter String path = GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory);

  impure function invertCase "Make one synthetic observation and append its inversion to the scene"
    input String path;
    input String caseName;
    input String title;
    input String description;
    input Real trueLens[7] "Radians";
    input Real trueSource[7] "Radians";
    input Boolean pointMass;
    input Real distances[3](each unit = "m");
    input Real apertureM(unit = "m");
    input Real fieldArcsec;
    input Real psfArcsec;
    input Real noiseFraction;
    input Integer pixels;
    input Integer seed;
    input Integer maxIterations;
    input Boolean writeLensRecord;
    input Real zLens;
    input Real zSource;
    input Real H0;
    input Real OmegaM;
    output Real ringRadius;
    output Real massInitialSolar;
    output Real thetaE;
    output Real thetaESigma;
    output Real massSolar;
    output Real massSigmaSolar;
    output Real trueThetaE;
    output Real trueMassSolar;
    output Real chi2Initial;
    output Real chi2Final;
    output Integer iterations;
    output Boolean converged;
  protected
    constant Real arcsec = GeneralRelativity.Constants.arcsec;
    Real u[pixels] = {(-fieldArcsec + (i - 0.5) * 2.0 * fieldArcsec / pixels) * arcsec for i in 1:pixels};
    Real noiseless[pixels, pixels];
    Real observed[pixels, pixels];
    Real noiseSigma;
    Real lens[7];
    Real sources[1, 7];
    Real lensSigma[7];
    Real massInApertureSolar;
    Real explained;
    Real correlation;
  algorithm
    noiseless := GeneralRelativity.Lensing.Inversion.renderImage(u, u, trueLens, {trueSource}, pointMass, 2, psfArcsec * arcsec);
    noiseSigma := noiseFraction * max(noiseless);
    observed := GeneralRelativity.Lensing.Inversion.syntheticImage(u, u, trueLens, {trueSource}, pointMass, noiseSigma, seed, 2, psfArcsec * arcsec);
    (ringRadius, massInitialSolar, lens, sources, lensSigma, thetaE, thetaESigma, massSolar, massSigmaSolar, massInApertureSolar, trueThetaE, trueMassSolar,
      chi2Initial, chi2Final, explained, correlation, iterations, converged) := GeneralRelativity.Visualization.writeLensInversion(path, caseName, title, description,
      u, u, observed, noiseSigma, pointMass, psfArcsec * arcsec, distances, apertureM, trueLens, trueSource, maxIterations, 4, 0.2,
      writeLensRecord, zLens, zSource, H0, OmegaM);
  end invertCase;

  parameter Real Dl[nCases](each fixed = false, each unit = "m");
  parameter Real Ds[nCases](each fixed = false, each unit = "m");
  parameter Real Dls[nCases](each fixed = false, each unit = "m");
  parameter Real trueB[nCases](each fixed = false, each unit = "rad") "True deflection scales";
  parameter Real ringRadius[nCases](each fixed = false);
  parameter Real thetaE[nCases](each fixed = false);
  parameter Real thetaESigma[nCases](each fixed = false);
  parameter Real trueThetaE[nCases](each fixed = false);
  parameter Real chi2Final[nCases](each fixed = false);
initial algorithm
  GeneralRelativity.Visualization.beginScene(path, "Lens inversion library", "cartesian", 1.0, "arcsec", 1.0);
  GeneralRelativity.Visualization.writeView(path, "lensfit");
  for k in 1:nCases loop
    Dl[k] := GeneralRelativity.Lensing.angularDiameterDistance(0.0, zLens[k], H0, OmegaM);
    Ds[k] := GeneralRelativity.Lensing.angularDiameterDistance(0.0, zSource[k], H0, OmegaM);
    Dls[k] := GeneralRelativity.Lensing.angularDiameterDistance(zLens[k], zSource[k], H0, OmegaM);
    trueB[k] := if pointMass[k] then GeneralRelativity.Lensing.einsteinRadiusPointMass(massSolar[k] * GeneralRelativity.Constants.GM_sun, Dl[k], Ds[k], Dls[k])
      else GeneralRelativity.Lensing.einsteinRadiusSIS(sqrt(massSolar[k] * GeneralRelativity.Constants.GM_sun / (GeneralRelativity.Constants.pi * apertureM)), Ds[k], Dls[k]);
    (ringRadius[k], massInitialSolar[k], thetaE[k], thetaESigma[k], massFitSolar[k], massSigmaSolar[k], trueThetaE[k], trueMassInRingSolar[k],
      chi2Initial[k], chi2Final[k], iterations[k], converged[k]) := invertCase(path, caseNames[k], titles[k], descriptions[k],
      {trueB[k], lensCentreArcsec[k, 1] * arcsec, lensCentreArcsec[k, 2] * arcsec, axisRatio[k], positionAngleDeg[k] * deg, shear[k], shearAngleDeg[k] * deg},
      {sourceArcsec[k, 1] * arcsec, sourceArcsec[k, 2] * arcsec, 1.0, sourceRadiusArcsec[k] * arcsec, sourceAxisRatio[k], sourceAngleDeg[k] * deg, sersicIndex[k]},
      pointMass[k], {Dl[k], Ds[k], Dls[k]}, apertureM, fieldArcsec[k], psfArcsec[k], noiseFraction[k], pixels, 13 * k + 5, maxIterations,
      k == 1, zLens[k], zSource[k], H0, OmegaM);
    trueThetaEArcsec[k] := trueThetaE[k] / arcsec;
    ringRadiusArcsec[k] := ringRadius[k] / arcsec;
    thetaEArcsec[k] := thetaE[k] / arcsec;
    thetaESigmaArcsec[k] := thetaESigma[k] / arcsec;
    massErrorSigmas[k] := (massFitSolar[k] - trueMassInRingSolar[k]) / max(massSigmaSolar[k], 1.0e-300);
    chi2PerPixel[k] := chi2Final[k] / pixels^2;
  end for;
equation
  for k in 1:nCases loop
    assert(not checkRecovery or chi2Final[k] < chi2Initial[k], "Iterating did not improve on the first simulation in case " + caseNames[k]);
    assert(not checkRecovery or (chi2PerPixel[k] > 0.7 and chi2PerPixel[k] < 1.5), "The fitted simulation does not match the image to within the noise in case " + caseNames[k]);
    assert(not checkRecovery or abs(thetaEArcsec[k] / trueThetaEArcsec[k] - 1.0) < 0.05, "The fitted Einstein radius is not within 5% of the truth in case " + caseNames[k]);
    assert(not checkRecovery or abs(massErrorSigmas[k]) < 4.0 or abs(massFitSolar[k] / trueMassInRingSolar[k] - 1.0) < 0.03,
      "The fitted mass is neither within 4 standard errors nor within 3% of the truth in case " + caseNames[k]);
  end for;
  annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 1.0e-8, Interval = 0.1),
    Documentation(info="<html>
<p>Eight synthetic observations that look different on the sky, each run through
<a href=\"modelica://GeneralRelativity.Visualization.writeLensInversion\">Visualization.writeLensInversion</a> (the steps of
<a href=\"modelica://GeneralRelativity.Lensing.Inversion\">Lensing.Inversion</a>) exactly as in
<a href=\"modelica://GeneralRelativity.Examples.LensInversion\">LensInversion</a>:</p>
<table>
<tr><th>Case</th><th>What the image shows</th></tr>
<tr><td>ring</td><td>Source almost on the axis of a round galaxy: a nearly complete Einstein ring.</td></tr>
<tr><td>double</td><td>Source well off the axis: a bright image and a faint counter-image.</td></tr>
<tr><td>arc</td><td>Flattened, sheared galaxy with the source near a caustic cusp: one long merged arc and a counter-image.</td></tr>
<tr><td>cross</td><td>Compact source near the centre of an elongated lens: four images (Einstein cross).</td></tr>
<tr><td>point</td><td>Point mass: two images on opposite sides, the inner one fainter and mirror-reversed.</td></tr>
<tr><td>noisy</td><td>Twice the seeing and four times the noise: same method, larger error bars.</td></tr>
<tr><td>group</td><td>A group-scale lens at z = 0.35 with a z = 3 source: arcs three arcseconds out.</td></tr>
<tr><td>small</td><td>A low-mass galaxy with an Einstein radius of a quarter of an arcsecond.</td></tr>
</table>
<p>Every case appends its observed, first-simulation, fitted and residual images and its fit summary to one scene, tagged with
the case name; the viewer's <b>lensfit</b> style shows them as a gallery to pick from. Per case the model reports the ring-radius
estimate (<code>massInitialSolar</code>), the fitted mass inside the critical curve with its standard error (<code>massFitSolar</code>,
<code>massSigmaSolar</code>), the truth (<code>trueMassInRingSolar</code>) and the error in standard errors
(<code>massErrorSigmas</code>). The assertions require every fit to reach the noise level, the Einstein radius within 5% and the
mass within four standard errors or 3%.</p>
<p>The images live inside a function, so the model holds only the per-case results and translates quickly. Add a case by
extending the parameter arrays and <code>nCases</code>; add a real image to a scene with
<a href=\"modelica://GeneralRelativity.Visualization.writeLensInversion\">writeLensInversion</a> as
<a href=\"modelica://GeneralRelativity.Examples.LensInversionFromFile\">LensInversionFromFile</a> does.</p>
</html>"));
end LensInversionLibrary;
