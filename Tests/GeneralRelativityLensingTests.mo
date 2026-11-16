within Tests;
model GeneralRelativityLensingTests "Gravitational lensing: cosmological distances, Einstein radii, point-mass images and the isothermal-ellipse deflection"
  function run
    output Boolean passed;
  protected
    constant Real Mpc = GeneralRelativity.Constants.Mpc;
    constant Real arcsec = GeneralRelativity.Constants.arcsec;
    Real Dl;
    Real Ds;
    Real Dls;
    Real thetaE;
    Real tp;
    Real tm;
    Real mp;
    Real mm;
    Real u;
    Real a[2];
    Real a2[2];
    Real ap[2];
    Real am[2];
    Real th[2];
    Real kappa;
    Real div;
    Real curl;
    Real h;
    Real b;
    Real q;
    Real sigma;
  algorithm
    // Distances in a flat universe, H0 = 70 km/s/Mpc, OmegaM = 0.3 (reference values from an independent
    // 200000-interval Simpson integration).
    Dl := GeneralRelativity.Lensing.angularDiameterDistance(0.0, 0.5);
    Ds := GeneralRelativity.Lensing.angularDiameterDistance(0.0, 2.0);
    Dls := GeneralRelativity.Lensing.angularDiameterDistance(0.5, 2.0);
    assert(abs(GeneralRelativity.Lensing.angularDiameterDistance(0.0, 1.0) / Mpc / 1651.914403 - 1.0) < 1.0e-8, "D_A(z = 1) is not 1651.91 Mpc");
    assert(abs(Dl / Mpc / 1259.083597 - 1.0) < 1.0e-8 and abs(Ds / Mpc / 1726.620691 - 1.0) < 1.0e-8, "D_A(0.5) or D_A(2) is wrong");
    assert(abs(Dls / Mpc / 1097.078893 - 1.0) < 1.0e-8, "D_A(0.5, 2) is wrong");
    // Nearby: Hubble's law D = c z / H0.
    assert(abs(GeneralRelativity.Lensing.angularDiameterDistance(0.0, 1.0e-4) / (GeneralRelativity.Constants.c * 1.0e-4 / (70.0e3 / Mpc)) - 1.0) < 1.0e-3,
      "Nearby distance is not cz/H0");
    // Curved universes (Hogg's formula; independent references) and the Einstein-de Sitter closed form
    // D_A(0, z) = 2 (c/H0) (1 - 1/sqrt(1 + z)) / (1 + z).
    assert(abs(GeneralRelativity.Lensing.angularDiameterDistance(0.0, 1.0, 70.0, 0.3, 0.6) / Mpc / 1619.958627 - 1.0) < 1.0e-8
      and abs(GeneralRelativity.Lensing.angularDiameterDistance(0.5, 2.0, 70.0, 0.3, 0.6) / Mpc / 1059.057545 - 1.0) < 1.0e-8,
      "Open-universe distances are wrong");
    assert(abs(GeneralRelativity.Lensing.angularDiameterDistance(0.0, 1.0, 70.0, 0.3, 0.8) / Mpc / 1686.527186 - 1.0) < 1.0e-8,
      "Closed-universe distance is wrong");
    assert(abs(GeneralRelativity.Lensing.angularDiameterDistance(0.0, 1.0, 70.0, 1.0, 0.0) / Mpc / 1254.388257 - 1.0) < 1.0e-8,
      "Einstein-de Sitter distance is not the closed form");
    assert(abs(GeneralRelativity.Lensing.angularDiameterDistance(0.5, 2.0, 70.0, 0.3, 0.7) - Dls) < 1.0e-6 * Dls, "Explicit flat case differs");

    // Einstein radii.
    thetaE := GeneralRelativity.Lensing.einsteinRadiusPointMass(1.0e12 * GeneralRelativity.Constants.GM_sun, Dl, Ds, Dls);
    assert(abs(thetaE / arcsec / 2.0272539 - 1.0) < 1.0e-6, "Point-mass Einstein radius of 1e12 solar masses is not 2.027 arcsec");
    sigma := sqrt(1.0e12 * GeneralRelativity.Constants.GM_sun / (GeneralRelativity.Constants.pi * 10.0e-3 * Mpc));
    assert(abs(sigma / 1000.0 / 370.00331 - 1.0) < 1.0e-6, "Velocity dispersion of 1e12 solar masses within 10 kpc is not 370 km/s");
    assert(abs(GeneralRelativity.Lensing.einsteinRadiusSIS(sigma, Ds, Dls) / arcsec / 2.5086826 - 1.0) < 1.0e-6, "SIS Einstein radius is not 2.509 arcsec");
    // The mass inside the Einstein ring is pi (Dl thetaE)^2 Sigma_cr, which for a point mass is the mass itself.
    assert(abs(GeneralRelativity.Constants.pi * (Dl * thetaE)^2 * GeneralRelativity.Lensing.criticalDensity(Dl, Ds, Dls) * GeneralRelativity.Constants.G
      / (1.0e12 * GeneralRelativity.Constants.GM_sun) - 1.0) < 1.0e-9, "Mass inside the Einstein ring is not the point mass");

    // Point-mass images solve the lens equation; magnifications obey the classic identities.
    u := 0.37;
    (tp, tm, mp, mm) := GeneralRelativity.Lensing.pointMassImages(u * thetaE, thetaE);
    assert(abs(tp - thetaE^2 / tp - u * thetaE) < 1.0e-12 * thetaE and abs(tm - thetaE^2 / tm - u * thetaE) < 1.0e-12 * thetaE,
      "Point-mass images do not solve the lens equation");
    assert(abs(mp + mm - 1.0) < 1.0e-12 and abs(abs(mp) + abs(mm) - (u^2 + 2.0) / (u * sqrt(u^2 + 4.0))) < 1.0e-12,
      "Point-mass magnifications are wrong");
    a := GeneralRelativity.Lensing.deflectionPointMass({tp, 0.0}, thetaE);
    assert(abs(a[1] - thetaE^2 / tp) < 1.0e-15 and abs(a[2]) < 1.0e-30, "Point-mass deflection is wrong");

    // Isothermal ellipse: at q = 1 it is the isothermal sphere (deflection b toward the centre everywhere);
    // in general half its divergence is the convergence q b / (2 psi) and its curl vanishes.
    b := 1.3 * arcsec;
    a := GeneralRelativity.Lensing.deflectionSIE({0.6 * arcsec, -0.8 * arcsec}, b, 1.0);
    assert(abs(a[1] - 0.6 * b) < 1.0e-15 and abs(a[2] + 0.8 * b) < 1.0e-15, "SIE with q = 1 is not the isothermal sphere");
    a2 := GeneralRelativity.Lensing.deflectionSIE({0.6 * arcsec, -0.8 * arcsec}, b, 0.9999999);
    assert(abs(a2[1] - a[1]) < 1.0e-6 * b and abs(a2[2] - a[2]) < 1.0e-6 * b, "SIE is not continuous at q = 1");
    q := 0.6;
    th := {0.7 * arcsec, 0.3 * arcsec};
    h := 1.0e-5 * arcsec;
    ap := GeneralRelativity.Lensing.deflectionSIE(th + {h, 0.0}, b, q);
    am := GeneralRelativity.Lensing.deflectionSIE(th - {h, 0.0}, b, q);
    div := (ap[1] - am[1]) / (2.0 * h);
    curl := (ap[2] - am[2]) / (2.0 * h);
    ap := GeneralRelativity.Lensing.deflectionSIE(th + {0.0, h}, b, q);
    am := GeneralRelativity.Lensing.deflectionSIE(th - {0.0, h}, b, q);
    div := div + (ap[2] - am[2]) / (2.0 * h);
    curl := curl - (ap[1] - am[1]) / (2.0 * h);
    kappa := GeneralRelativity.Lensing.convergenceSIE(th, b, q);
    assert(abs(0.5 * div / kappa - 1.0) < 1.0e-6, "Half the divergence of the SIE deflection is not its convergence");
    assert(abs(curl) < 1.0e-6 * abs(div), "SIE deflection is not a gradient (non-zero curl)");

    // External shear: the deflection of the potential (g1 (x^2 - y^2) + 2 g2 x y)/2 is linear and traceless.
    a := GeneralRelativity.Lensing.deflectionShear({1.0, 0.0}, 0.1, 0.0);
    a2 := GeneralRelativity.Lensing.deflectionShear({0.0, 1.0}, 0.1, 0.0);
    assert(abs(a[1] - 0.1) < 1.0e-15 and abs(a2[2] + 0.1) < 1.0e-15 and abs(a[2]) < 1.0e-15, "Shear deflection is wrong");
    a := GeneralRelativity.Lensing.deflectionShear({1.0, 0.0}, 0.1, GeneralRelativity.Constants.pi / 4.0);
    assert(abs(a[1]) < 1.0e-15 and abs(a[2] - 0.1) < 1.0e-15, "Shear position angle is wrong");
    passed := true;
  end run;

  parameter Boolean passed = run();
equation
  assert(passed, "GeneralRelativity lensing assertions did not complete");
  annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 1.0e-8, Interval = 0.001));
end GeneralRelativityLensingTests;
