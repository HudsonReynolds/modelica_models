within Tests;
model GeneralRelativityCosmologyTests "FLRW cosmology: Friedmann solutions, distances, recombination, spherical Bessel functions"
  function run
    input GeneralRelativity.Cosmology.Parameters par "Planck 2018";
    output Boolean passed;
  protected
    import C = GeneralRelativity.Cosmology;
    C.Parameters mL(omegaBh2 = 0.0, omegaCh2 = 0.3 * 0.7^2, h = 0.7, TCMB = 1.0e-3, Neff = 0.0) "Matter + Lambda (radiation negligible)";
    C.Parameters rad(h = 0.7, omegaBh2 = 1.0e-12, omegaCh2 = 1.0e-12, TCMB = 30.0, Neff = 3.046) "Radiation-dominated early on";
    Real j[2501];
    Real j5[6];
    Real j60[61];
    Real x;
    Real t;
    Real a;
    Real H0;
    Real Og, On, Ob, Oc, OL;
    Real z[6000], aa[6000], eta[6000], xe[6000], dtau[6000], tau[6000], g[6000], kD[6000];
    Real zStar, etaStar, zReio, eta0, rs, sumG, chi;
  algorithm
    // Spherical Bessel functions: closed forms and reference values (scipy.special.spherical_jn).
    x := 2.3;
    j := C.sphericalBessel(2500, x);
    assert(abs(j[1] - sin(x) / x) < 1.0e-14 and abs(j[2] - (sin(x) / x^2 - cos(x) / x)) < 1.0e-14
      and abs(j[3] - ((3.0 / x^2 - 1.0) * sin(x) / x - 3.0 * cos(x) / x^2)) < 1.0e-14, "j_0, j_1, j_2 are not the closed forms");
    j5 := C.sphericalBessel(5, 3.7);
    assert(abs(j5[6] - 0.03861365693381352) < 1.0e-14, "j_5(3.7) is wrong");
    j := C.sphericalBessel(2500, 1100.25);
    assert(abs(j[1001] + 0.0011509281162058798) < 1.0e-12, "j_1000(1100.25) is wrong");
    j := C.sphericalBessel(2500, 2000.0);
    assert(abs(j[2001] - 0.0009586857291362243) < 1.0e-12, "j_2000(2000) is wrong (transition region)");
    j60 := C.sphericalBessel(60, 10.0);
    assert(abs(j60[51] / 2.2306960232186405e-31 - 1.0) < 1.0e-10, "j_50(10) is wrong (far below the turning point)");
    // Recurrence j_{l-1} + j_{l+1} = (2l+1)/x j_l.
    j := C.sphericalBessel(2500, 2700.0);
    assert(abs(j[2501] - 0.0005816901537739074) < 1.0e-12, "j_2500(2700) is wrong");
    assert(abs(j[1500] + j[1502] - 3001.0 / 2700.0 * j[1501]) < 1.0e-13, "The recurrence does not hold");

    // Flat matter + Lambda: t(a) = 2 / (3 H0 sqrt(OL)) asinh(sqrt(OL/Om) a^1.5).
    (H0, Og, On, Ob, Oc, OL) := C.densities(mL);
    for a1 in {0.1, 0.5, 1.0, 2.0} loop
      t := C.cosmicTime(a1, mL);
      assert(abs(t / (2.0 / (3.0 * H0 * sqrt(OL)) * C.MpcGyr * Modelica.Math.asinh(sqrt(OL / (Ob + Oc)) * a1^1.5)) - 1.0) < 1.0e-6,
        "cosmicTime is not the matter + Lambda closed form");
      assert(abs(C.scaleFactorAtTime(t, mL) / a1 - 1.0) < 1.0e-9, "scaleFactorAtTime does not invert cosmicTime");
    end for;
    // The same universe: comoving distance and Lensing.angularDiameterDistance (which has no radiation).
    chi := C.comovingDistance(1.0 / 1.5, 1.0, mL);
    assert(abs(chi / 1.5 / (GeneralRelativity.Lensing.angularDiameterDistance(0.0, 0.5, 70.0, 0.3, OL) / GeneralRelativity.Constants.Mpc) - 1.0) < 1.0e-5,
      "comovingDistance disagrees with Lensing.angularDiameterDistance");
    // Radiation era: a = sqrt(2 H0 sqrt(Or) t), conformal time eta = a / (H0 sqrt(Or)).
    (H0, Og, On, Ob, Oc, OL) := C.densities(rad);
    a := 1.0e-6;
    assert(abs(C.cosmicTime(a, rad) / (a^2 / (2.0 * H0 * sqrt(Og + On)) * C.MpcGyr) - 1.0) < 1.0e-6, "The radiation era is not a ~ t^1/2");
    assert(abs(C.conformalTime(a, rad) / (a / (H0 * sqrt(Og + On))) - 1.0) < 1.0e-6, "The radiation era is not a ~ eta");
    // Event horizon of pure de Sitter would be c/H; with matter it is larger, and finite.
    assert(C.eventHorizon(1.0, par) > GeneralRelativity.Constants.c / 1000.0 / C.hubble(1.0, par) and C.eventHorizon(1.0, par) < 6000.0,
      "The event horizon today is not between c/H0 and 6 Gpc");

    // Planck 2018 cosmology: age, last scattering, the acoustic scale, the damping scale.
    assert(abs(C.cosmicTime(1.0, par) - 13.79) < 0.03, "The age is not 13.79 +- 0.03 Gyr");
    (z, aa, eta, xe, dtau, tau, g, kD, zStar, etaStar, zReio) := C.thermalHistory(par);
    eta0 := C.conformalTime(1.0, par);
    rs := C.soundHorizon(1.0 / (1.0 + zStar), par);
    assert(abs(zStar - 1089.9) < 10.0, "Last scattering is not at z* = 1090");
    assert(abs(zReio - 7.67) < 0.3, "The reionization redshift for tau = 0.054 is not 7.7");
    assert(abs(rs / 144.43 - 1.0) < 0.005, "The sound horizon at last scattering is not 144.4 Mpc");
    assert(abs(100.0 * rs / (eta0 - etaStar) / 1.04110 - 1.0) < 0.005, "100 theta* is not 1.0411");
    // The visibility function integrates to 1 (every CMB photon last scattered somewhere).
    sumG := 0.0;
    for i in 2:6000 loop
      sumG := sumG + 0.5 * (g[i] + g[i - 1]) * (eta[i] - eta[i - 1]);
    end for;
    assert(abs(sumG - (1.0 - exp(-tau[1]))) < 1.0e-3, "The visibility function does not integrate to 1");
    // Residual ionization after recombination ~ 2e-4, and full ionization before it.
    assert(xe[5990] > 1.0 and xe[1] > 1.07, "Reionized / early universe is not fully ionized");
    passed := true;
  end run;

  parameter GeneralRelativity.Cosmology.Parameters par;
  // Evaluated in compiled code when the simulation starts (as a parameter binding, the frontend would
  // interpret the recombination integration, which takes over a minute).
  discrete Boolean passed(start = false, fixed = true);
algorithm
  when initial() then
    passed := run(par);
    assert(passed, "GeneralRelativity cosmology assertions did not complete");
  end when;
  annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 1.0e-8, Interval = 0.5));
end GeneralRelativityCosmologyTests;
