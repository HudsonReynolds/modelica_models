within GeneralRelativity.Examples;
model CosmicMicrowaveBackground "Acoustic peaks of the cosmic microwave background from linear perturbations on the FLRW metric, compared with Planck 2018; writes the viewer's CMB scene"
  parameter GeneralRelativity.Cosmology.Parameters par "Lambda-CDM parameters (Planck 2018 by default)";
  parameter String sceneFile = "cmb.ndjson" "Scene file name for Visualization/threejs/index.html";
  parameter String sceneDirectory = "" "Output folder; empty uses GR_SCENE_DIR, else ~/Documents/Visualization, else the working directory";
  discrete Real summary[24] "See Cosmology.writeCmbScene";
  discrete Real zStar "Redshift of last scattering (peak of the visibility function)";
  discrete Real thetaStar "100 theta*: angular size of the sound horizon at last scattering";
  discrete Real rsStar(unit = "Mpc") "Comoving sound horizon at last scattering";
  discrete Real lPeak1 "Multipole of the first acoustic peak (model)";
  discrete Real DPeak1(unit = "uK2") "Height of the first peak (model)";
  discrete Real lPeak1Data "Multipole of the first peak in the Planck data";
  discrete Real DPeak1Data(unit = "uK2") "Height of the first peak in the Planck data";
  discrete Real chi2 "Chi-square against the Planck spectrum (diagonal errors; indicative)";
  discrete Real rmsVsBestFit "RMS fractional difference from Planck's best-fit spectrum, l = 30..1500";
algorithm
  when initial() then
    summary := GeneralRelativity.Cosmology.writeCmbScene(par, GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory));
    zStar := summary[1];
    rsStar := summary[4];
    thetaStar := summary[5];
    lPeak1 := summary[9];
    DPeak1 := summary[10];
    chi2 := summary[17];
    rmsVsBestFit := summary[19];
    lPeak1Data := summary[20];
    DPeak1Data := summary[21];
    assert(abs(zStar - 1089.9) < 10.0, "Last scattering is not at z = 1090");
    assert(abs(thetaStar / 1.04110 - 1.0) < 0.005, "The acoustic angular scale differs from Planck's 100 theta* = 1.0411 by more than 0.5%");
    assert(abs(lPeak1 / lPeak1Data - 1.0) < 0.03, "The first acoustic peak is not where Planck measures it");
    assert(abs(DPeak1 / DPeak1Data - 1.0) < 0.15, "The first acoustic peak height differs from Planck's by more than 15%");
  end when;
  annotation(experiment(StartTime = 0.0, StopTime = 1.0, Interval = 0.5),
    Documentation(info="<html>
<p>The temperature anisotropies of the cosmic microwave background, computed from first principles on the FLRW
metric with the Planck 2018 parameters and compared with the Planck 2018 measurements. Everything is computed
when the simulation starts (a few seconds):</p>
<ol>
<li><b>Recombination</b> (<a href=\"modelica://GeneralRelativity.Cosmology.thermalHistory\">thermalHistory</a>):
the free-electron fraction from Peebles' atom, the visibility function (last scattering at
<code>z<sub>*</sub> = 1090</code>) and the diffusion-damping scale.</li>
<li><b>Sound waves</b> (<a href=\"modelica://GeneralRelativity.Cosmology.acousticSources\">acousticSources</a>):
450 Fourier modes of the photon-baryon plasma, dark matter, neutrinos and the metric potentials, each evolved
from the radiation era. Modes that have entered the horizon oscillate: the plasma's pressure pushes back against
gravity, with sound speed <code>c/&radic;(3(1 + R))</code>.</li>
<li><b>Projection onto the sky</b> (<a href=\"modelica://GeneralRelativity.Cosmology.cmbSpectrum\">cmbSpectrum</a>):
the line-of-sight integral with spherical Bessel functions turns the state of each mode at last scattering into
multipoles; the phase of each mode when the CMB was released sets the acoustic peaks.</li>
</ol>
<p>Results (Planck 2018 parameters, massless neutrinos): <code>z<sub>*</sub> = 1089.8</code>, sound horizon
<code>r<sub>s</sub> = 144.4 Mpc</code>, <code>100&theta;<sub>*</sub> = 1.0396</code> (Planck: 1.0411), damping scale
<code>k<sub>D</sub> = 0.141/Mpc</code>, drag epoch <code>z = 1061</code>. Acoustic peaks at <code>l = 223, 538, 816</code>
against <code>223, 535, 819</code> in the data; heights within 4%; the spectrum is within 2.5% (rms, <code>l = 30..1500</code>)
of Planck's best-fit spectrum, the rest being mostly CMB lensing, which smooths Planck's peaks and is not modelled.
The chi-square against the 111 data points with diagonal errors is about 700 (the data's per-cent errors at
high <code>l</code> resolve that lensing).</p>
<p>The scene <code>cmb.ndjson</code> opens in the viewer's <b>cmb</b> style: the power spectrum with the data and its
contributions, the recombination history, and the transfer functions of a few hundred modes at 40 times from
the radiation era to after recombination, from which the viewer draws the expanding sound shell of a point-like
perturbation (Eisenstein, Seo &amp; White 2007), frozen at the sound horizon when the baryons are released.</p>
</html>"));
end CosmicMicrowaveBackground;
