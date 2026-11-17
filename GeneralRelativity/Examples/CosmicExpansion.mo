within GeneralRelativity.Examples;
model CosmicExpansion "Time evolution of a Lambda-CDM universe on the FLRW metric: the Friedmann equation integrated in cosmic time, horizons, galaxies and light rays; writes the viewer's expanding-universe scene"
  import Cosmo = GeneralRelativity.Cosmology;
  parameter Cosmo.Parameters par "Lambda-CDM parameters (Planck 2018 by default)";
  parameter String sceneFile = "cosmic_expansion.ndjson" "Scene file name for Visualization/threejs/index.html";
  parameter String sceneDirectory = "" "Output folder; empty uses GR_SCENE_DIR, else ~/Documents/Visualization, else the working directory";
  parameter Real galaxyRedshifts[4] = {0.5, 1.0, 2.0, 5.0} "Galaxies whose light reaches us today from these redshifts";
  parameter Real galaxyAngles[4] = {20.0, 110.0, 200.0, 290.0} "Their directions in the plane, degrees";
  parameter Real sampleInterval(unit = "Gyr") = 0.05 "Time between scene samples";
  final parameter Real t0(unit = "Gyr") = 1.0e-6 "Start (the experiment's StartTime): deep in the radiation era";
  final parameter Real H0Gyr(unit = "1/Gyr") = Cosmo.hubble(1.0, par) * 1000.0 / GeneralRelativity.Constants.Mpc * Cosmo.Gyr "Hubble constant in 1/Gyr";
  // Quantities from the quadratures, fixed at the start.
  discrete Real age(unit = "Gyr") "Age of the universe: a = 1";
  discrete Real eta0(unit = "Gly") "Conformal time today: the comoving radius of the observable universe";
  discrete Real zStar "Redshift of last scattering";
  discrete Real tStar(unit = "Gyr") "Cosmic time of last scattering";
  discrete Real etaStar(unit = "Gly") "Conformal time of last scattering";
  discrete Real chiGalaxy[4](each unit = "Gly") "Comoving distances of the galaxies";
  discrete Real epochs[12](each start = 0.0, each fixed = true) "See Cosmology.writeExpansionRecords";
  discrete Real aCheck(start = 0.0, fixed = true) "a(t) from the quadrature at the last check";
  Real lnA(start = log(Cosmo.scaleFactorAtTime(1.0e-6, par)), fixed = true) "ln a, integrated from the Friedmann equation";
  Real eta(start = Cosmo.conformalTime(Cosmo.scaleFactorAtTime(1.0e-6, par), par) * Cosmo.MpcGyr, fixed = true) "Conformal time, Gly: d eta/dt = c/a";
  Real a = exp(lnA) "Scale factor";
  Real z = 1.0 / a - 1.0 "Redshift of light emitted now and received today";
  Real H(unit = "km/(s.Mpc)") = Cosmo.hubble(a, par) "Hubble rate";
  Real T(unit = "K") = par.TCMB / a "Temperature of the CMB";
  Real q = Cosmo.decelerationParameter(a, par) "Deceleration parameter (negative: accelerating)";
  Real hubbleRadius(unit = "Gly") = GeneralRelativity.Constants.c / 1000.0 / H * Cosmo.MpcGyr "Proper Hubble radius c/H";
  Real particleHorizon(unit = "Gly") = a * eta "Proper distance to the particle horizon";
  Real dGalaxy[4](each unit = "Gly") = a * chiGalaxy "Proper distances of the galaxies";
  Real dCmbPhoton(unit = "Gly") = a * (if time >= tStar then eta0 - eta else eta0 - etaStar)
    "Proper distance from us of the CMB photon we receive today (before last scattering: of the plasma that will emit it)";
  Real dPhotonToday(unit = "Gly") = if time >= age then a * (eta - eta0) else 0.0 "Proper distance of a photon we send out today";
  GeneralRelativity.Visualization.SceneRecorder scene(
    spherical = false,
    metricName = "",
    sceneFile = sceneFile,
    sceneDirectory = sceneDirectory,
    title = "Expanding universe (Planck 2018 LCDM)",
    referenceRadius = 1.0,
    lengthUnit = "Gly",
    c = 1.0,
    rOuter = 60.0,
    writeEmbedding = false,
    writeProfile = false,
    writeSlice = false,
    nTracks = 6,
    trackNames = cat(1, {"galaxy z = " + String(galaxyRedshifts[i]) for i in 1:4}, {"CMB photon", "photon sent today"}),
    trackColors = {"#ffb454", "#7ee0a1", "#c792ea", "#ff6b8a", "#5fb3ff", "#f9e2af"},
    trackDescriptions = cat(1, {"A comoving galaxy whose light, emitted at z = " + String(galaxyRedshifts[i]) + ", reaches us today" for i in 1:4},
      {"The CMB photon we receive today: emitted at last scattering, its proper distance first grows (space expands faster than it approaches), then shrinks",
       "A photon we send out today: it will never get beyond the event horizon"}),
    nScalars = 2,
    scalarNames = fill({"z", "distance"}, 6),
    sampleInterval = sampleInterval,
    style = "cosmos",
    trackX = cat(1, {{time, dGalaxy[i] * cos(galaxyAngles[i] * GeneralRelativity.Constants.pi / 180.0), dGalaxy[i] * sin(galaxyAngles[i] * GeneralRelativity.Constants.pi / 180.0), 0.0} for i in 1:4},
      {{time, dCmbPhoton, 0.0, 0.0}, {time, -dPhotonToday * cos(0.3), -dPhotonToday * sin(0.3), 0.0}}),
    trackLambda = fill(time, 6),
    trackScalars = cat(1, {{z, dGalaxy[i]} for i in 1:4}, {{z, dCmbPhoton}, {z, dPhotonToday}}));
equation
  // The Friedmann equation, d(ln a)/dt = H(a), and conformal time.
  der(lnA) = H * 1000.0 / GeneralRelativity.Constants.Mpc * Cosmo.Gyr;
  der(eta) = 1.0 / a;
algorithm
  when initial() then
    age := Cosmo.cosmicTime(1.0, par);
    eta0 := Cosmo.conformalTime(1.0, par) * Cosmo.MpcGyr;
    zStar := Cosmo.lastScatteringRedshift(par);
    tStar := Cosmo.cosmicTime(1.0 / (1.0 + zStar), par);
    etaStar := Cosmo.conformalTime(1.0 / (1.0 + zStar), par) * Cosmo.MpcGyr;
    for i in 1:4 loop
      chiGalaxy[i] := eta0 - Cosmo.conformalTime(1.0 / (1.0 + galaxyRedshifts[i]), par) * Cosmo.MpcGyr;
    end for;
  end when;
  // The records are written once the scene file has been started (by the recorder's initial algorithm).
  when time >= t0 + 1.0e-4 then
    epochs := Cosmo.writeExpansionRecords(par, GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory), zStar);
  end when;
  // Check the integrated a(t) against the quadrature t(a), every 2 Gyr and at a = 1.
  when sample(1.0, 2.0) then
    aCheck := Cosmo.scaleFactorAtTime(time, par);
    assert(abs(a / aCheck - 1.0) < 1.0e-6, "The integrated scale factor disagrees with the quadrature of the Friedmann equation");
  end when;
  when a >= 1.0 then
    assert(abs(time / age - 1.0) < 1.0e-6, "a = 1 is not reached at the age of the universe");
  end when;
  annotation(experiment(StartTime = 1.0e-6, StopTime = 60.0, Interval = 0.01, Tolerance = 1.0e-10),
    Documentation(info="<html>
<p>The flat FLRW metric <code>ds&sup2; = &minus;c&sup2;dt&sup2; + a(t)&sup2;dx&sup2;</code> with the Planck 2018 contents:
the Friedmann equation <code>d(ln a)/dt = H(a)</code> is integrated from a thousand years after the big bang
(radiation era) through matter domination into the accelerating, Lambda-dominated future at 60 Gyr, together
with conformal time <code>d&eta;/dt = c/a</code>. Every 2 Gyr <code>a(t)</code> is checked against the quadrature
<a href=\"modelica://GeneralRelativity.Cosmology.cosmicTime\">cosmicTime</a>, and <code>a = 1</code> must be
reached at the age (13.81 Gyr).</p>
<p>Epochs (in <code>epochs</code>): matter-radiation equality at <code>z = 3402</code> (51,000 years), last scattering
at <code>z = 1090</code> (372,000 years), matter-Lambda equality at <code>z = 0.30</code> (10.3 Gyr), accelerating
expansion since <code>z = 0.63</code> (7.7 Gyr). Today the Hubble radius is 14.5 Gly, the particle horizon (the
edge of the observable universe) 46.2 Gly, the event horizon 16.7 Gly, and the CMB was emitted by matter now 45.3
Gly away.</p>
<p>Distances are in Gly and times in Gyr (<code>c = 1 Gly/Gyr</code>). The scene <code>cosmic_expansion.ndjson</code>
opens in the viewer's <b>cosmos</b> style. Its tracks are proper distances from us: four galaxies whose light,
emitted at <code>z = 0.5, 1, 2, 5</code>, reaches us today (they recede with the expansion); the CMB photon we
receive today (before last scattering, the plasma that will emit it), whose distance first grows because
space expands faster than it approaches, then shrinks to zero today; and a photon we send out today, which
never gets beyond the event horizon. The <code>cosmology</code> and <code>expansion</code> records hold the
epochs and the expansion history (a, z, H, T, horizons, deceleration parameter) for the viewer's diagrams.</p>
</html>"));
end CosmicExpansion;
