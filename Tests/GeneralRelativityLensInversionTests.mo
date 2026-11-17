within Tests;
model GeneralRelativityLensInversionTests "Lens inversion: mass from the Einstein radius, Sersic sources, the image simulation, the radius read off an image, the residual metrics, the mass inside an elliptical critical curve, and the fit with and without a prior and from several starts"
  function run
    output Boolean passed;
  protected
    constant Real arcsec = GeneralRelativity.Constants.arcsec;
    constant Real pi = GeneralRelativity.Constants.pi;
    constant Integer n = 33 "Pixels per axis; odd so that one pixel centre is at the origin";
    Real Dl;
    Real Ds;
    Real Dls;
    Real GM;
    Real thetaE;
    Real sigma;
    Real u[n];
    Real lens[7];
    Real src[1, 7];
    Real src2[1, 7];
    Real image[n, n];
    Real image2[n, n];
    Real observed[n, n];
    Real a[2];
    Real a2[2];
    Real bn;
    Real r;
    Real c[2];
    Real spread;
    Real ss;
    Real chi2;
    Real explained;
    Real corr;
    Real noise;
    Real sumD;
    Real sumD2;
    Real lens0[7];
    Real src0[1, 7];
    Real fit[7];
    Real fitSrc[1, 7];
    Real history[21];
    Real lensHistory[21, 7];
    Integer iterations;
    Boolean converged;
    Real fitSigma[7];
    Real fitSrcSigma[1, 7];
    Real prior[7];
    Integer mid = div(n + 1, 2);
    Real twoSources[2, 7];
    Real reference;
    Real theta[2];
    Real beta[2];
    Real thetaEff;
    Real thetaEffSigma;
    Real GMfit;
    Real GMfitSigma;
    Real cov[7, 7];
    Real kappaSum;
    Real q;
    Real b;
    Real psi;
    Real h;
    Real ringRadius;
    Real ringCentre[2];
    Real lensStart[7];
    Real srcStart[1, 7];
    Real history2[4 + 30 + 1];
    Real lensHistory2[4 + 30 + 1, 7];
    Real fitCov[7, 7];
    constant Integer nk = 401 "Grid for the convergence integral";
  algorithm
    Dl := GeneralRelativity.Lensing.angularDiameterDistance(0.0, 0.5);
    Ds := GeneralRelativity.Lensing.angularDiameterDistance(0.0, 2.0);
    Dls := GeneralRelativity.Lensing.angularDiameterDistance(0.5, 2.0);
    // Mass inside the Einstein radius inverts the point-mass Einstein radius; the SIS dispersion inverts its radius.
    GM := 1.0e12 * GeneralRelativity.Constants.GM_sun;
    thetaE := GeneralRelativity.Lensing.einsteinRadiusPointMass(GM, Dl, Ds, Dls);
    assert(abs(GeneralRelativity.Lensing.massInEinsteinRadius(thetaE, Dl, Ds, Dls) / GM - 1.0) < 1.0e-12, "massInEinsteinRadius does not invert einsteinRadiusPointMass");
    sigma := 250.0e3;
    assert(abs(GeneralRelativity.Lensing.velocityDispersionSIS(GeneralRelativity.Lensing.einsteinRadiusSIS(sigma, Ds, Dls), Ds, Dls) / sigma - 1.0) < 1.0e-12,
      "velocityDispersionSIS does not invert einsteinRadiusSIS");

    // Sersic source: I(re) = Ie, I(0) = Ie exp(bn); q > 1 is the same ellipse with the axes swapped.
    src := [0.3 * arcsec, -0.2 * arcsec, 2.0, 0.25 * arcsec, 0.6, 0.4, 1.0];
    bn := 2.0 - 1.0 / 3.0 + 4.0 / 405.0 + 46.0 / 25515.0;
    assert(abs(GeneralRelativity.Lensing.Inversion.sersicSource(src[1, 1:2] + 0.25 * arcsec * {cos(0.4), sin(0.4)}, src[1, :]) - 2.0) < 1.0e-12, "Sersic brightness at re is not Ie");
    assert(abs(GeneralRelativity.Lensing.Inversion.sersicSource(src[1, 1:2], src[1, :]) - 2.0 * exp(bn)) < 1.0e-9, "Sersic central brightness is not Ie exp(bn)");
    src2 := src;
    src2[1, 5] := 1.0 / 0.6;
    src2[1, 6] := 0.4 + 0.5 * pi;
    for k in 1:5 loop
      a := {0.17 * k - 0.4, 0.31 * k - 0.9} * arcsec;
      assert(abs(GeneralRelativity.Lensing.Inversion.sersicSource(a, src[1, :]) - GeneralRelativity.Lensing.Inversion.sersicSource(a, src2[1, :])) < 1.0e-12,
        "A source with q > 1 is not the swapped ellipse");
    end for;

    // Deflection of the parameter vector: an isothermal sphere deflects by thetaE everywhere; a point mass as deflectionPointMass.
    lens := {1.2 * arcsec, 0.1 * arcsec, -0.05 * arcsec, 1.0, 0.3, 0.0, 0.0};
    a := GeneralRelativity.Lensing.Inversion.deflectionAt({0.9 * arcsec, 0.4 * arcsec}, lens);
    assert(abs(sqrt(a * a) - 1.2 * arcsec) < 1.0e-15, "SIS deflection of the parameter vector is not thetaE");
    a := GeneralRelativity.Lensing.Inversion.deflectionAt({0.9 * arcsec, 0.4 * arcsec}, lens, true);
    a2 := GeneralRelativity.Lensing.deflectionPointMass({0.8 * arcsec, 0.45 * arcsec}, 1.2 * arcsec);
    assert(abs(a[1] - a2[1]) < 1.0e-18 and abs(a[2] - a2[2]) < 1.0e-18, "Point-mass deflection of the parameter vector is wrong");

    // Rendering: no lens shows the source itself; a point mass with a centred round source makes a symmetric ring.
    u := {(i - mid) * 0.1 * arcsec for i in 1:n};
    src := [0.0, 0.0, 1.0, 0.2 * arcsec, 1.0, 0.0, 1.0];
    image := GeneralRelativity.Lensing.Inversion.renderImage(u, u, {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0}, src, false, 1);
    assert(abs(image[mid, mid] - exp(bn)) < 1.0e-9, "Without a lens the image is not the source");
    assert(abs(image[mid + 2, mid] - exp(-bn * (0.2 * arcsec / (0.2 * arcsec) - 1.0))) < 1.0e-9, "Without a lens the image profile is not the Sersic profile");
    image := GeneralRelativity.Lensing.Inversion.renderImage(u, u, {1.2 * arcsec, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0}, src, true, 2);
    for i in 1:n loop
      for j in 1:n loop
        assert(abs(image[i, j] - image[n + 1 - i, j]) < 1.0e-12 and abs(image[i, j] - image[i, n + 1 - j]) < 1.0e-12 and abs(image[i, j] - image[j, i]) < 1.0e-12,
          "A centred source behind a point mass does not make a symmetric ring");
      end for;
    end for;
    // The ring is brightest near thetaE and empty at the centre (the source's light is spread into the ring).
    assert(image[mid + 12, mid] > 10.0 * image[mid + 2, mid] and image[mid + 12, mid] > 10.0 * image[mid + 16, mid], "The point-mass ring is not at the Einstein radius");
    // The point-spread function conserves flux and blurs.
    image2 := GeneralRelativity.Lensing.Inversion.renderImage(u, u, {1.2 * arcsec, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0}, src, true, 2, 0.1 * arcsec);
    assert(abs(sum(image2) / sum(image) - 1.0) < 0.02, "The point-spread function does not conserve flux");
    assert(max(image2) < max(image), "The point-spread function does not blur");

    // Synthetic observation: zero noise reproduces the render; the noise has the requested scatter.
    observed := GeneralRelativity.Lensing.Inversion.syntheticImage(u, u, {1.2 * arcsec, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0}, src, true, 0.0, 3, 2);
    assert(max(abs(observed - image)) < 1.0e-15, "syntheticImage without noise differs from renderImage");
    noise := 0.05 * max(image);
    observed := GeneralRelativity.Lensing.Inversion.syntheticImage(u, u, {1.2 * arcsec, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0}, src, true, noise, 3, 2);
    sumD := sum(observed - image) / n^2;
    sumD2 := sum((observed - image) .^ 2) / n^2;
    assert(abs(sumD) < 0.15 * noise and abs(sqrt(sumD2 - sumD^2) / noise - 1.0) < 0.1, "The synthetic noise does not have the requested scatter");

    // The radius read off the image: exact for an isothermal sphere even with an offset source (the two images lie on a circle of radius thetaE).
    (thetaE, c, spread) := GeneralRelativity.Lensing.Inversion.einsteinRadiusFromImage(u, u, image, 0.2);
    assert(abs(thetaE / (1.2 * arcsec) - 1.0) < 0.03 and sqrt(c * c) < 0.05 * arcsec, "The circle fitted to a point-mass ring has the wrong radius or centre");
    src := [0.15 * arcsec, 0.1 * arcsec, 1.0, 0.1 * arcsec, 1.0, 0.0, 1.0];
    lens := {1.2 * arcsec, 0.1 * arcsec, -0.05 * arcsec, 1.0, 0.0, 0.0, 0.0};
    image := GeneralRelativity.Lensing.Inversion.renderImage(u, u, lens, src, false, 2);
    (thetaE, c, spread) := GeneralRelativity.Lensing.Inversion.einsteinRadiusFromImage(u, u, image, 0.2);
    assert(abs(thetaE / (1.2 * arcsec) - 1.0) < 0.05, "The circle fitted to the two images of an isothermal sphere is not the Einstein radius");
    assert(sqrt((c[1] - 0.15 * arcsec)^2 + (c[2] - 0.1 * arcsec)^2) < 0.08 * arcsec, "The circle fitted to the two images is not centred on the source position (lens centre plus source offset)");
    assert(spread < 0.15 * thetaE, "The two images of an isothermal sphere do not lie near one circle");
    // The centroid of a symmetric ring is its centre.
    image := GeneralRelativity.Lensing.Inversion.renderImage(u, u, {1.2 * arcsec, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0}, [0.0, 0.0, 1.0, 0.2 * arcsec, 1.0, 0.0, 1.0], true, 2);
    c := GeneralRelativity.Lensing.Inversion.brightnessCentroid(u, u, image, 0.2);
    assert(sqrt(c * c) < 1.0e-12 * arcsec, "The centroid of a symmetric ring is not its centre");

    // Residual metrics.
    (ss, chi2, explained, corr) := GeneralRelativity.Lensing.Inversion.imageResidual(image, image, 0.1);
    assert(ss == 0.0 and chi2 == 0.0 and abs(explained - 1.0) < 1.0e-15 and abs(corr - 1.0) < 1.0e-12, "A perfect match does not score perfectly");
    (ss, chi2, explained, corr) := GeneralRelativity.Lensing.Inversion.imageResidual(image, zeros(n, n), 0.1);
    assert(abs(ss - sum(image .^ 2)) < 1.0e-12 * ss and abs(chi2 - ss / 0.01) < 1.0e-9 * chi2 and abs(explained) < 1.0e-15 and corr == 0.0, "An empty simulation does not score zero");
    (ss, chi2, explained, corr) := GeneralRelativity.Lensing.Inversion.imageResidual(image, 3.0 * image + fill(0.5, n, n), 0.1);
    assert(abs(corr - 1.0) < 1.0e-12, "Correlation is not invariant to scale and offset");

    // The fit: an isothermal ellipse with an offset exponential source, 2% noise, started from the circle fit.
    lens := {1.2 * arcsec, 0.1 * arcsec, -0.05 * arcsec, 0.8, 20.0 * pi / 180.0, 0.0, 0.0};
    src := [0.15 * arcsec, 0.1 * arcsec, 1.0, 0.12 * arcsec, 1.0, 0.0, 1.0];
    image := GeneralRelativity.Lensing.Inversion.renderImage(u, u, lens, src, false, 2);
    noise := 0.02 * max(image);
    observed := GeneralRelativity.Lensing.Inversion.syntheticImage(u, u, lens, src, false, noise, 11, 2);
    (thetaE, c, spread) := GeneralRelativity.Lensing.Inversion.einsteinRadiusFromImage(u, u, observed, 0.2);
    lens0 := {thetaE, c[1], c[2], 0.9, 0.0, 0.0, 0.0};
    src0 := [c[1], c[2], max(observed), 0.3 * thetaE, 1.0, 0.0, 1.0];
    (fit, fitSrc, history, lensHistory, iterations, converged, fitSigma, fitSrcSigma) := GeneralRelativity.Lensing.Inversion.fitLens(u, u, observed, noise, lens0, src0,
      {true, true, true, true, true, false, false}, [true, true, true, true, false, false, false], maxIterations = 20);
    assert(converged, "The lens fit did not converge in 20 iterations");
    for k in 2:21 loop
      assert(history[k] <= history[k - 1], "The fit cost increased");
    end for;
    assert(history[21] < history[1] and history[21] > 0.7 * n^2 and history[21] < 1.5 * n^2, "The fit does not reach the noise level");
    assert(abs(fit[1] / (1.2 * arcsec) - 1.0) < 0.01, "The fitted Einstein radius is not within 1% of the truth");
    assert(abs(fit[1] - 1.2 * arcsec) < 3.0 * fitSigma[1] + 0.005 * arcsec, "The fitted Einstein radius is not within its standard error of the truth");
    assert(abs(fit[4] - 0.8) < 0.05 and abs(fit[5] - 20.0 * pi / 180.0) < 0.1, "The fitted lens shape is wrong");
    assert(sqrt((fit[2] - 0.1 * arcsec)^2 + (fit[3] + 0.05 * arcsec)^2) < 0.03 * arcsec, "The fitted lens position is wrong");
    assert(sqrt((fitSrc[1, 1] - 0.15 * arcsec)^2 + (fitSrc[1, 2] - 0.1 * arcsec)^2) < 0.03 * arcsec, "The fitted source position is wrong");
    assert(fitSigma[6] == 0.0 and fitSigma[7] == 0.0 and fitSrcSigma[1, 5] == 0.0, "Fixed parameters report a standard error");
    assert(abs(lensHistory[21, 1] - fit[1]) < 1.0e-15 and abs(lensHistory[1, 1] - thetaE) < 1.0e-15, "The lens history does not start at the guess and end at the fit");
    assert(abs(fit[5] - lensHistory[21, 5]) < 1.0e-12 or abs(abs(fit[5] - lensHistory[21, 5]) - pi) < 1.0e-9, "The canonical position angle differs from the history by other than pi");
    // A tight prior (the Kalman update with outside knowledge) holds the Einstein radius where the prior says.
    prior := zeros(7);
    prior[1] := 1.0e-4 * arcsec;
    lens0[1] := 1.1 * 1.2 * arcsec;
    (fit, fitSrc, history, lensHistory, iterations, converged, fitSigma, fitSrcSigma) := GeneralRelativity.Lensing.Inversion.fitLens(u, u, observed, noise, lens0, src0,
      {true, true, true, true, true, false, false}, [true, true, true, true, false, false, false], prior, maxIterations = 20);
    assert(abs(fit[1] / (1.1 * 1.2 * arcsec) - 1.0) < 0.01, "A tight prior does not hold the Einstein radius");
    assert(fitSigma[1] < 1.1e-4 * arcsec, "The standard error under a tight prior is not the prior's");

    // The fast renderer computes exactly the reference functions: deflectionAt and sersicSource, per pixel, with an ellipse lens,
    // shear, and two sources (one with q > 1).
    lens := {1.1 * arcsec, 0.07 * arcsec, -0.04 * arcsec, 0.7, 0.6, 0.06, 1.1};
    twoSources := [0.12 * arcsec, 0.05 * arcsec, 1.0, 0.15 * arcsec, 0.6, 0.4, 1.0; -0.2 * arcsec, 0.1 * arcsec, 0.5, 0.1 * arcsec, 1.6, 1.3, 2.5];
    image := GeneralRelativity.Lensing.Inversion.renderImage(u, u, lens, twoSources, false, 1);
    for i in 1:n loop
      for j in 1:n loop
        theta := {u[i], u[j]};
        beta := theta - GeneralRelativity.Lensing.Inversion.deflectionAt(theta, lens);
        reference := GeneralRelativity.Lensing.Inversion.sersicSource(beta, twoSources[1, :]) + GeneralRelativity.Lensing.Inversion.sersicSource(beta, twoSources[2, :]);
        assert(abs(image[i, j] - reference) < 1.0e-10 * max(reference, 1.0e-3), "renderImage differs from deflectionAt and sersicSource");
      end for;
    end for;
    image := GeneralRelativity.Lensing.Inversion.renderImage(u, u, lens, twoSources, true, 1);
    theta := {u[5], u[20]};
    beta := theta - GeneralRelativity.Lensing.Inversion.deflectionAt(theta, lens, true);
    assert(abs(image[5, 20] - GeneralRelativity.Lensing.Inversion.sersicSource(beta, twoSources[1, :]) - GeneralRelativity.Lensing.Inversion.sersicSource(beta, twoSources[2, :]))
      < 1.0e-10, "renderImage differs from the references for a point mass");

    // Mass inside the critical curve of an isothermal ellipse: the convergence integrated over psi < q b is pi q b^2 (mean 1),
    // and einsteinMass reports it through the equal-area radius b sqrt(q).
    q := 0.6;
    b := 1.0;
    h := 2.0 * b / q / (nk - 1);
    kappaSum := 0.0;
    for i in 1:nk loop
      for j in 1:nk loop
        theta := {-b / q + (i - 1) * h + 0.5 * h * 0.0, -b / q + (j - 1) * h};
        psi := sqrt(q^2 * theta[1]^2 + theta[2]^2);
        // Skip the grid point at the centre by distance, not by psi > 0: with fused multiply-adds (Apple silicon) it lands
        // about 1e-16 off the origin, and its 1/psi would swamp the sum.
        if psi < q * b and psi > 0.5 * q * h then
          kappaSum := kappaSum + GeneralRelativity.Lensing.convergenceSIE(theta, b, q) * h^2;
        end if;
      end for;
    end for;
    assert(abs(kappaSum / (GeneralRelativity.Constants.pi * q * b^2) - 1.0) < 0.02, "The convergence inside the critical ellipse of an SIE is not pi q b^2");
    lens := {1.2 * arcsec, 0.0, 0.0, 0.64, 0.3, 0.05, 0.2};
    cov := zeros(7, 7);
    (thetaEff, thetaEffSigma, GMfit, GMfitSigma) := GeneralRelativity.Lensing.Inversion.einsteinMass(lens, cov, false, Dl, Ds, Dls);
    assert(abs(thetaEff / (0.8 * 1.2 * arcsec) - 1.0) < 1.0e-12 and thetaEffSigma == 0.0 and GMfitSigma == 0.0, "The equal-area Einstein radius is not b sqrt(q)");
    assert(abs(GMfit / GeneralRelativity.Lensing.massInEinsteinRadius(1.2 * arcsec, Dl, Ds, Dls) - 0.64) < 1.0e-12, "The mass inside the critical ellipse is not q times that of the sphere with the same b");
    (thetaEff, thetaEffSigma, GMfit, GMfitSigma) := GeneralRelativity.Lensing.Inversion.einsteinMass(lens, cov, true, Dl, Ds, Dls);
    assert(abs(thetaEff / (1.2 * arcsec) - 1.0) < 1.0e-12, "A point mass's Einstein radius is not b");
    // Error propagation: d(thetaE)/thetaE = db/b + dq/(2q); the mass error is twice the relative radius error, and correlation counts.
    cov[1, 1] := (0.01 * 1.2 * arcsec)^2;
    cov[4, 4] := 0.02^2;
    (thetaEff, thetaEffSigma, GMfit, GMfitSigma) := GeneralRelativity.Lensing.Inversion.einsteinMass(lens, cov, false, Dl, Ds, Dls);
    assert(abs(thetaEffSigma / thetaEff - sqrt(0.01^2 + (0.02 / (2.0 * 0.64))^2)) < 1.0e-12 and abs(GMfitSigma / GMfit - 2.0 * thetaEffSigma / thetaEff) < 1.0e-12,
      "The Einstein-radius and mass errors are not propagated from b and q");
    cov[1, 4] := -0.9 * sqrt(cov[1, 1] * cov[4, 4]);
    cov[4, 1] := cov[1, 4];
    (thetaEff, thetaEffSigma, GMfit, GMfitSigma) := GeneralRelativity.Lensing.Inversion.einsteinMass(lens, cov, false, Dl, Ds, Dls);
    assert(thetaEffSigma / thetaEff < sqrt(0.01^2 + (0.02 / (2.0 * 0.64))^2), "Anti-correlated b and q do not shrink the Einstein-radius error");

    // invert on a double: one bright image, one faint counter-image, where a single fit started at the circle through the images can converge to
    // the wrong lens (that circle is centred on the source); the multi-start must recover it.
    lens := {1.0 * arcsec, -0.1 * arcsec, 0.05 * arcsec, 0.85, 70.0 * pi / 180.0, 0.0, 0.0};
    src := [0.45 * arcsec, 0.2 * arcsec, 1.0, 0.1 * arcsec, 0.7, 2.0, 1.0];
    image := GeneralRelativity.Lensing.Inversion.renderImage(u, u, lens, src, false, 2, 0.05 * arcsec);
    noise := 0.02 * max(image);
    observed := GeneralRelativity.Lensing.Inversion.syntheticImage(u, u, lens, src, false, noise, 23, 2, 0.05 * arcsec);
    (ringRadius, ringCentre, lensStart, srcStart, fit, fitSrc, history2, lensHistory2, iterations, converged, fitSigma, fitSrcSigma, fitCov) :=
      GeneralRelativity.Lensing.Inversion.invert(u, u, observed, noise, false, 0.05 * arcsec, 30, 4);
    assert(history2[1] > history2[4 + 30 + 1], "invert did not improve on its start");
    assert(history2[4 + 30 + 1] < 1.5 * n^2, "invert did not reach the noise level on a double");
    (thetaEff, thetaEffSigma, GMfit, GMfitSigma) := GeneralRelativity.Lensing.Inversion.einsteinMass(fit, fitCov, false, Dl, Ds, Dls);
    assert(abs(thetaEff / (sqrt(0.85) * arcsec) - 1.0) < 0.05, "invert did not recover the Einstein radius of a double within 5%");
    assert(sqrt((fit[2] + 0.1 * arcsec)^2 + (fit[3] - 0.05 * arcsec)^2) < 0.1 * arcsec, "invert did not find the lens position of a double");
    passed := true;
  end run;

  parameter Boolean passed = run();
equation
  assert(passed, "GeneralRelativity lens-inversion assertions did not complete");
  annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 1.0e-8, Interval = 0.001));
end GeneralRelativityLensInversionTests;
