within GeneralRelativity.Examples;
model GalaxyLensing "Gravitational lensing by a foreground galaxy: cosmological distances, Einstein radii and a microlensing light curve, with a lensed-sky scene for the viewer"
  parameter Real massSolar = 5.0e11 "Lens mass in solar masses: projected within the aperture (galaxy) or total (point)";
  parameter Real apertureKpc = 10.0 "Aperture radius for the galaxy mass";
  parameter Real zLens = 0.5 "Lens redshift";
  parameter Real zSource = 2.0 "Source redshift";
  parameter Real H0 = 70.0 "Hubble constant, km/s/Mpc";
  parameter Real OmegaM = 0.3 "Matter density (flat universe)";
  parameter Boolean pointMass = false "true: point-mass lens; false: isothermal-ellipse galaxy with external shear";
  parameter Real axisRatio = 0.75 "Galaxy axis ratio";
  parameter Real positionAngleDeg = 30.0 "Galaxy position angle";
  parameter Real shear = 0.05 "External shear from the galaxy's environment";
  parameter Real shearAngleDeg = 110.0 "Shear position angle";
  parameter Real u0 = 0.2 "Closest approach of the moving source, in Einstein radii (microlensing light curve)";
  parameter String sceneFile = "galaxy_lensing.ndjson" "Scene file name";
  parameter String sceneDirectory = "" "Output folder; empty uses GR_SCENE_DIR, else ~/Documents/Visualization";
  final parameter Real Dl(unit = "m") = GeneralRelativity.Lensing.angularDiameterDistance(0.0, zLens, H0, OmegaM) "Distance to the lens";
  final parameter Real Ds(unit = "m") = GeneralRelativity.Lensing.angularDiameterDistance(0.0, zSource, H0, OmegaM) "Distance to the source";
  final parameter Real Dls(unit = "m") = GeneralRelativity.Lensing.angularDiameterDistance(zLens, zSource, H0, OmegaM) "Lens-to-source distance";
  final parameter Real GM(unit = "m3/s2") = massSolar * GeneralRelativity.Constants.GM_sun;
  final parameter Real sigma(unit = "m/s") = sqrt(GM / (GeneralRelativity.Constants.pi * apertureKpc * 1.0e-3 * GeneralRelativity.Constants.Mpc))
    "Velocity dispersion of the isothermal galaxy with that mass inside the aperture: M(<R) = pi sigma^2 R / G";
  final parameter Real thetaEPoint(unit = "rad") = GeneralRelativity.Lensing.einsteinRadiusPointMass(GM, Dl, Ds, Dls);
  final parameter Real thetaEGalaxy(unit = "rad") = GeneralRelativity.Lensing.einsteinRadiusSIS(sigma, Ds, Dls);
  final parameter Real thetaE(unit = "rad") = if pointMass then thetaEPoint else thetaEGalaxy "Einstein radius of the chosen model";
  final parameter Real SigmaCr(unit = "kg/m2") = GeneralRelativity.Lensing.criticalDensity(Dl, Ds, Dls);
  final parameter Real massInEinsteinRadiusSolar = GeneralRelativity.Constants.pi * (Dl * thetaE)^2 * SigmaCr * GeneralRelativity.Constants.G / GeneralRelativity.Constants.GM_sun
    "Projected mass inside the Einstein radius, pi (Dl thetaE)^2 Sigma_cr, in solar masses";
  // Microlensing: a source drifting behind a point mass. Modelica time is in units of the Einstein-radius crossing time.
  Real u "Source offset in Einstein radii";
  Real thetaPlus(unit = "rad");
  Real thetaMinus(unit = "rad");
  Real muPlus;
  Real muMinus;
  Real totalMagnification "|mu+| + |mu-|, the microlensing light curve";
equation
  u = sqrt(u0^2 + (time - 0.5)^2 / 0.15^2);
  (thetaPlus, thetaMinus, muPlus, muMinus) = GeneralRelativity.Lensing.pointMassImages(u * thetaEPoint, thetaEPoint);
  totalMagnification = abs(muPlus) + abs(muMinus);
  assert(abs(u * thetaEPoint - (thetaPlus - thetaEPoint^2 / thetaPlus)) < 1.0e-12 * thetaEPoint
    and abs(u * thetaEPoint - (thetaMinus - thetaEPoint^2 / thetaMinus)) < 1.0e-12 * thetaEPoint,
    "Point-mass images do not solve the lens equation");
  assert(abs(totalMagnification - (u^2 + 2.0) / (u * sqrt(u^2 + 4.0))) < 1.0e-9 * totalMagnification,
    "Total magnification is not (u^2 + 2)/(u sqrt(u^2 + 4))");
  assert(abs(muPlus + muMinus - 1.0) < 1.0e-9 * totalMagnification, "Signed magnifications do not sum to 1");
  assert(abs(massInEinsteinRadiusSolar / massSolar - (if pointMass then 1.0 else thetaE * Dl / (apertureKpc * 1.0e-3 * GeneralRelativity.Constants.Mpc))) < 1.0e-9,
    "Mass inside the Einstein radius is inconsistent with the lens model");
initial algorithm
  GeneralRelativity.Visualization.beginScene(GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory), "Galaxy lensing", "cartesian", 1.0, "arcsec", 1.0);
  GeneralRelativity.Visualization.writeView(GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory), "lensing");
  GeneralRelativity.Visualization.writeLens(GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory),
    if pointMass then "point" else "galaxy", massSolar, apertureKpc, zLens, zSource, H0, OmegaM, axisRatio, positionAngleDeg,
    shear, shearAngleDeg, thetaE / GeneralRelativity.Constants.arcsec, {Dl, Ds, Dls} / GeneralRelativity.Constants.Mpc);
  annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 1.0e-8, Interval = 0.002),
    Documentation(info="<html>
<p>A galaxy at redshift <code>z<sub>l</sub></code> lenses galaxies at <code>z<sub>s</sub></code>. The model
computes the angular-diameter distances of a flat &Lambda;CDM universe, the critical density, the
velocity dispersion of an isothermal galaxy with the given mass inside the aperture
(<code>M = &pi;&sigma;&sup2;R/G</code>), and the Einstein radius for an isothermal galaxy or a point mass of the
same mass. With the defaults (<code>5&times;10<sup>11</sup> M<sub>&#9737;</sub></code> within 10 kpc,
<code>z = 0.5 &rarr; 2</code>) &sigma; = 262 km/s and &theta;<sub>E</sub> = 1.25 arcsec, a typical galaxy-scale
strong lens.</p>
<p>As a simulation it follows a source drifting behind a point mass (time in units of the
Einstein-crossing time): both image positions, their magnifications and the microlensing light curve,
asserting the lens equation, <code>|&mu;<sub>+</sub>| + |&mu;<sub>&minus;</sub>| = (u&sup2; + 2)/(u&radic;(u&sup2; + 4))</code>
and <code>&mu;<sub>+</sub> + &mu;<sub>&minus;</sub> = 1</code>.</p>
<p>The scene opens in the viewer's <b>lensing</b> style: a telescope-like image of background galaxies
seen through the lens, with sliders for the lens mass, shape, shear and redshifts; critical curves,
caustics in a source-plane inset, and convergence or magnification maps. The viewer recomputes the
distances live with the same formulas and shows the Modelica Einstein radius for comparison.</p>
</html>"));
end GalaxyLensing;
