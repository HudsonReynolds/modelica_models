within Tests;
model GeneralRelativityMetricTests "Exact metrics, electromagnetic and polytropic sources, PPN tests, gauge invariance and spherical harmonics"
  function kerrNewmanPotential "Kerr-Newman vector potential A_a = -(Q r / Sigma) (1, 0, 0, -a sin^2 theta), geometrized Gaussian"
    input Real x[4] "{ct, r, theta, phi}";
    input Real a "Spin length";
    input Real Q "Charge length";
    output Real A[4] "Covariant potential";
  protected
    Real Sigma = x[2]^2 + a^2 * cos(x[3])^2;
  algorithm
    A := -Q * x[2] / Sigma * {1.0, 0.0, 0.0, -a * sin(x[3])^2};
  end kerrNewmanPotential;

  function gaugeCoefficient "Cubic coefficient C[b, i, j, k] of the test gauge vector xi_b = C x x x + q x x + l x"
    input Integer b;
    input Integer i;
    input Integer j;
    input Integer k;
    output Real C;
  algorithm
    C := 0.1 * sin(b + 2.0 * i + 3.0 * j + 5.0 * k);
  end gaugeCoefficient;

  function gaugeGradient "Analytic d_a xi_b of the cubic test gauge vector"
    input Real x[4];
    output Real dxi[4, 4] "dxi[a, b] = d_a xi_b";
  algorithm
    for a in 1:4 loop
      for b in 1:4 loop
        dxi[a, b] := sum((gaugeCoefficient(b, a, j, k) + gaugeCoefficient(b, j, a, k) + gaugeCoefficient(b, j, k, a)) * x[j] * x[k]
          for j in 1:4, k in 1:4)
          + sum((cos(b + a * j) + cos(b + j * a)) * x[j] for j in 1:4) + 0.3 * a - 0.2 * b;
      end for;
    end for;
  end gaugeGradient;

  function gaugeThirdDerivative "Analytic d_k d_l d_a xi_b of the cubic test gauge vector"
    output Real dddxi[4, 4, 4, 4] "dddxi[k, l, a, b]";
  algorithm
    for k in 1:4 loop
      for l in 1:4 loop
        for a in 1:4 loop
          for b in 1:4 loop
            dddxi[k, l, a, b] := gaugeCoefficient(b, k, l, a) + gaugeCoefficient(b, k, a, l) + gaugeCoefficient(b, l, k, a)
              + gaugeCoefficient(b, l, a, k) + gaugeCoefficient(b, a, k, l) + gaugeCoefficient(b, a, l, k);
          end for;
        end for;
      end for;
    end for;
  end gaugeThirdDerivative;

  function pureGaugeSecondDifference "Central second differences of the pure-gauge h_ab = -(d_a xi_b + d_b xi_a)"
    input Real x[4];
    input Real s "Step";
    output Real ddh[4, 4, 4, 4];
  protected
    Real e[4, 4] = identity(4);
    Real hpp[4, 4];
    Real hpm[4, 4];
    Real hmp[4, 4];
    Real hmm[4, 4];
    Real dd[4, 4];
  algorithm
    for k in 1:4 loop
      for l in 1:4 loop
        hpp := GeneralRelativity.Perturbation.gaugeTransform(zeros(4, 4), gaugeGradient(x + s * e[k, :] + s * e[l, :]));
        hpm := GeneralRelativity.Perturbation.gaugeTransform(zeros(4, 4), gaugeGradient(x + s * e[k, :] - s * e[l, :]));
        hmp := GeneralRelativity.Perturbation.gaugeTransform(zeros(4, 4), gaugeGradient(x - s * e[k, :] + s * e[l, :]));
        hmm := GeneralRelativity.Perturbation.gaugeTransform(zeros(4, 4), gaugeGradient(x - s * e[k, :] - s * e[l, :]));
        dd := (hpp - hpm - hmp + hmm) / (4.0 * s^2);
        for a in 1:4 loop
          for b in 1:4 loop
            ddh[k, l, a, b] := dd[a, b];
          end for;
        end for;
      end for;
    end for;
  end pureGaugeSecondDifference;

  function run
    output Boolean passed;
  protected
    constant Real pi = 3.141592653589793;
    constant Real step[4] = {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3};
    constant Real eta[4, 4] = [-1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
    constant Real arcsec = 180.0 * 3600.0 / 3.141592653589793 "Arcseconds per radian";
    GeneralRelativity.Tensors.Curvature k;
    GeneralRelativity.Tensors.Curvature k2;
    Real x[4];
    Real g[4, 4];
    Real g2[4, 4];
    Real F[4, 4];
    Real T[4, 4];
    Real A[4];
    Real Aplus[4];
    Real Aminus[4];
    Real dA[4, 4] "dA[a, b] = d_a A_b";
    Real xs[4];
    Real u[4];
    Real rs;
    Real rQ;
    Real spin;
    Real Lambda;
    Real r;
    Real rho;
    Real H;
    Real accel;
    Real scale;
    Real kappa;
    Real U;
    Real err00;
    Real errSpace;
    Real E[3];
    Real B[3];
    Real mu0;
    Real cLight;
    Real value;
    Real value2;
    Real dvalue;
    Real b;
    Real r1;
    Real r2;
    Real chord;
    GeneralRelativity.Perturbation.PPNParameters ppnBD;
    Real dddxi[4, 4, 4, 4];
    Real ddh[4, 4, 4, 4];
    Real ddhFD[4, 4, 4, 4];
    Real ddhNew[4, 4, 4, 4];
    Real R[4, 4, 4, 4];
    Real R2[4, 4, 4, 4];
    Real hessian[3, 3];
    Real Glin[4, 4];
    Real laplacian;
    Real h[4, 4];
    Real hNew[4, 4];
    Real omega[4, 4];
    Real waveK[4];
    Real phase;
    Real nodes[6];
    Real weights[6];
    Real nodes10[10];
    Real weights10[10];
    Real Yre[25, 72];
    Real Yim[25, 72];
    Real Yreal[25, 72];
    Real w[72];
    Real gram;
    Real gramReal;
    Real err;
    Real re;
    Real im;
    Real re2;
    Real im2;
    Real re3;
    Real im3;
    Real sumSq;
    Real quadRe;
    Real quadIm;
    Real theta;
    Real phi;
    Real cg;
    Real cg2;
    Real w1;
    Real w2;
    Integer idx;
    Integer idx2;
    Integer lmTable[25, 2];
    Integer gauntCases[4, 6];
  algorithm
    // ---------------------------------------------------------------- Reissner-Nordstrom
    rs := 2.0;
    rQ := 0.8;
    x := {0.0, 10.0, 1.1, 0.4};
    k := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.reissnerNordstrom, x, {rs, rQ}, step);
    assert(k.valid > 0.5, "Reissner-Nordstrom metric was not inverted");
    // Finite differences with step 1e-3 leave an O(h^2) error of about 1e-8 in R (the same as for Schwarzschild).
    assert(abs(k.ricciScalar) < 5.0e-8, "Reissner-Nordstrom Ricci scalar is not zero (electrovacuum is traceless): " + String(k.ricciScalar));
    // Coulomb field F_tr = Q/r^2 in geometrized Gaussian units, Q = rQ.
    F := zeros(4, 4);
    F[1, 2] := rQ / x[2]^2;
    F[2, 1] := -F[1, 2];
    T := GeneralRelativity.StressEnergy.electromagnetic(k.g, F);
    assert(Tests.Assertions.maxAbsMatrix(k.einstein - 8.0 * pi * T) < 3.0e-8,
      "Reissner-Nordstrom Einstein tensor is not 8 pi T of the Coulomb field: "
      + String(Tests.Assertions.maxAbsMatrix(k.einstein - 8.0 * pi * T)));
    assert(abs(T[1, 1] - (1.0 - rs / x[2] + rQ^2 / x[2]^2) * rQ^2 / (8.0 * pi * x[2]^4)) < 1.0e-14,
      "Coulomb T_tt is not f Q^2 / (8 pi r^4)");
    assert(abs(GeneralRelativity.StressEnergy.trace(k.gInv, T)) < 1.0e-15, "Electromagnetic stress-energy is not traceless");
    assert(Tests.Assertions.maxAbsMatrix(GeneralRelativity.Metrics.reissnerNordstrom(x, {rs, 0.0})
      - GeneralRelativity.Metrics.schwarzschild(x, {rs})) < 1.0e-15, "Reissner-Nordstrom with rQ = 0 is not Schwarzschild");

    // ---------------------------------------------------------------- Kerr-Newman
    spin := 0.6;
    rQ := 0.5;
    x := {0.0, 6.0, 1.1, 0.4};
    assert(Tests.Assertions.maxAbsMatrix(GeneralRelativity.Metrics.kerrNewman(x, {rs, 0.0, rQ})
      - GeneralRelativity.Metrics.reissnerNordstrom(x, {rs, rQ})) < 1.0e-15, "Kerr-Newman with a = 0 is not Reissner-Nordstrom");
    assert(Tests.Assertions.maxAbsMatrix(GeneralRelativity.Metrics.kerrNewman(x, {rs, spin, 0.0})
      - GeneralRelativity.Metrics.kerr(x, {rs, spin})) < 1.0e-15, "Kerr-Newman with rQ = 0 is not Kerr");
    k := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.kerrNewman, x, {rs, spin, rQ}, step);
    assert(abs(k.ricciScalar) < 1.0e-7, "Kerr-Newman Ricci scalar is not zero");
    // F = dA from the analytic potential by central differences (first derivatives, step 1e-4).
    for i in 1:4 loop
      xs := x;
      xs[i] := x[i] + 1.0e-4;
      Aplus := kerrNewmanPotential(xs, spin, rQ);
      xs[i] := x[i] - 1.0e-4;
      Aminus := kerrNewmanPotential(xs, spin, rQ);
      dA[i, :] := (Aplus - Aminus) / 2.0e-4;
    end for;
    F := dA - transpose(dA);
    T := GeneralRelativity.StressEnergy.electromagnetic(k.g, F);
    assert(Tests.Assertions.maxAbsMatrix(k.einstein - 8.0 * pi * T) < 1.0e-7,
      "Kerr-Newman Einstein tensor is not 8 pi T of F = dA");

    // ---------------------------------------------------------------- Schwarzschild-de Sitter
    Lambda := 0.01;
    x := {0.0, 6.0, 1.1, 0.4};
    k := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.schwarzschildDeSitter, x, {rs, Lambda}, step);
    // Compare the mixed tensor R^a_b = Lambda delta^a_b so the r^2-sized angular components do not dominate.
    assert(Tests.Assertions.maxAbsMatrix(k.gInv * k.ricci - Lambda * identity(4)) < 1.0e-7,
      "Schwarzschild-de Sitter: R_ab is not Lambda g_ab");
    assert(abs(k.ricciScalar - 4.0 * Lambda) < 1.0e-7, "Schwarzschild-de Sitter: R is not 4 Lambda");
    assert(Tests.Assertions.maxAbsMatrix(GeneralRelativity.StressEnergy.einsteinResidual(k.einstein, k.g, zeros(4, 4), 1.0, 1.0, Lambda))
      < 1.0e-7, "Schwarzschild-de Sitter does not solve G + Lambda g = 0");
    // Static de Sitter is maximally symmetric: K = 8 Lambda^2 / 3.
    k := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.schwarzschildDeSitter, x, {0.0, Lambda}, step);
    assert(abs(k.kretschmann / (8.0 * Lambda^2 / 3.0) - 1.0) < 1.0e-5, "de Sitter Kretschmann scalar is not 8 Lambda^2 / 3");

    // ---------------------------------------------------------------- FLRW
    // Dust, a = (ct/ct0)^(2/3): G_ab = 8 pi rho u_a u_b with rho = 3 H^2 / (8 pi), H = n / (ct).
    x := {1.0, 0.3, -0.2, 0.5};
    k := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.flrw, x, {0.7, 2.0 / 3.0}, step);
    H := (2.0 / 3.0) / x[1];
    u := {1.0, 0.0, 0.0, 0.0};
    T := GeneralRelativity.StressEnergy.dust(k.g, u, 3.0 * H^2 / (8.0 * pi));
    assert(abs(k.g[2, 2] - (1.0 / 0.7)^(4.0 / 3.0)) < 1.0e-14, "FLRW scale factor is not (ct/ct0)^n");
    assert(Tests.Assertions.maxAbsMatrix(k.einstein - 8.0 * pi * T) < 5.0e-7,
      "Einstein-de Sitter FLRW Einstein tensor is not 8 pi times dust with rho = 3 H^2 / (8 pi)");
    // Radiation, n = 1/2: traceless Ricci and a radiation fluid with the same Friedmann density.
    k := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.flrw, x, {0.7, 0.5}, step);
    H := 0.5 / x[1];
    assert(abs(k.ricciScalar) < 1.0e-7, "Radiation-dominated FLRW has nonzero Ricci scalar");
    T := GeneralRelativity.StressEnergy.radiation(k.g, u, 3.0 * H^2 / (8.0 * pi), 1.0);
    assert(Tests.Assertions.maxAbsMatrix(k.einstein - 8.0 * pi * T) < 1.0e-7,
      "Radiation FLRW Einstein tensor is not 8 pi times a radiation fluid");
    // de Sitter expansion a = exp(H (ct - ct0)): R_ab = 3 H^2 g_ab.
    H := 0.5;
    k := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.flrw, x, {0.7, 0.0, H}, step);
    assert(Tests.Assertions.maxAbsMatrix(k.gInv * k.ricci - 3.0 * H^2 * identity(4)) < 1.0e-6, "Exponential FLRW: R_ab is not 3 H^2 g_ab");
    // Curved FLRW: Friedmann constraint G_00 = 3 (H^2 + k/a^2) and G^r_r = G^theta_theta = -(2 a''/a + H^2 + k/a^2).
    x := {1.0, 0.9, 1.1, 0.4};
    k := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.flrwCurved, x, {0.8, 0.5, 0.3}, step);
    scale := (x[1] / 0.8)^0.5;
    H := 0.5 / x[1];
    accel := 0.5 * (0.5 - 1.0) / x[1]^2; // a''/a
    kappa := 0.3 / scale^2;
    assert(abs(k.einstein[1, 1] - 3.0 * (H^2 + kappa)) < 1.0e-7, "Curved FLRW violates the Friedmann constraint");
    assert(abs(k.gInv[2, 2] * k.einstein[2, 2] + (2.0 * accel + H^2 + kappa)) < 2.0e-6
      and abs(k.gInv[3, 3] * k.einstein[3, 3] + (2.0 * accel + H^2 + kappa)) < 2.0e-6,
      "Curved FLRW spatial Einstein tensor is not -(2 a''/a + H^2 + k/a^2) or not isotropic");

    // ---------------------------------------------------------------- Isotropic Schwarzschild
    x := {0.0, 3.0, 2.0, 1.0};
    rho := sqrt(14.0);
    r := rho * (1.0 + rs / (4.0 * rho))^2;
    k := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.isotropicSchwarzschild, x, {rs}, step);
    assert(Tests.Assertions.maxAbsMatrix(k.ricci) < 1.0e-7, "Isotropic Schwarzschild is not vacuum");
    assert(abs(k.kretschmann / (12.0 * rs^2 / r^6) - 1.0) < 1.0e-5,
      "Isotropic Schwarzschild Kretschmann is not 12 rs^2 / r^6 at the areal radius");
    assert(abs(k.g[1, 1] + (1.0 - rs / r)) < 1.0e-14, "Isotropic g_00 is not -(1 - rs/r) at the matching areal radius");

    // ---------------------------------------------------------------- PPN metric
    rs := 1.0;
    x := {0.0, 600.0, 0.0, 800.0}; // r = 1000 rs
    U := rs / 2000.0;
    g := GeneralRelativity.Metrics.ppn(x, {rs, 1.0, 1.0});
    g2 := GeneralRelativity.Metrics.isotropicSchwarzschild(x, {rs});
    err00 := g2[1, 1] - g[1, 1];
    errSpace := g2[2, 2] - g[2, 2];
    assert(abs(err00 / (1.5 * U^3) - 1.0) < 1.0e-2, "PPN (gamma = beta = 1) g_00 differs from Schwarzschild by more than 3 U^3 / 2");
    assert(abs(errSpace / (1.5 * U^2) - 1.0) < 1.0e-2, "PPN g_ij differs from Schwarzschild by more than 3 U^2 / 2");
    assert(abs(g[2, 3]) < 1.0e-15 and abs(g[3, 3] - g[2, 2]) < 1.0e-15 and abs(g[1, 2]) < 1.0e-15,
      "PPN spatial metric is not isotropic and static");
    g := GeneralRelativity.Metrics.ppn(x, {rs, 0.5, 2.0});
    assert(abs(g[3, 3] - (1.0 + 2.0 * 0.5 * U)) < 1.0e-15 and abs(g[1, 1] - (-1.0 + 2.0 * U - 2.0 * 2.0 * U^2)) < 1.0e-15,
      "PPN metric does not carry gamma and beta");
    ppnBD := GeneralRelativity.Perturbation.bransDicke(500.0);
    assert(abs(ppnBD.gamma - 501.0 / 502.0) < 1.0e-15 and ppnBD.beta == 1.0 and ppnBD.alpha1 == 0.0,
      "Brans-Dicke PPN parameters are wrong");
    assert(abs(GeneralRelativity.Perturbation.nordtvedtEta(GeneralRelativity.Perturbation.PPNParameters())) < 1.0e-15,
      "Nordtvedt parameter is not zero in general relativity");
    assert(abs(GeneralRelativity.Perturbation.nordtvedtEta(ppnBD) - 1.0 / 502.0) < 1.0e-15,
      "Brans-Dicke Nordtvedt parameter is not 1/(2 + omega)");

    // ---------------------------------------------------------------- Classical tests
    value := GeneralRelativity.Perturbation.ppnPerihelionAdvance(GeneralRelativity.Constants.GM_sun, 5.7909e10, 0.2056)
      * (36525.0 / 87.969) * arcsec;
    assert(abs(value / 42.98 - 1.0) < 5.0e-3, "Mercury's perihelion advance is not 42.98 arcsec per century");
    value2 := GeneralRelativity.Perturbation.ppnPerihelionAdvance(GeneralRelativity.Constants.GM_sun, 5.7909e10, 0.2056, 0.5, 1.0)
      * (36525.0 / 87.969) * arcsec;
    assert(abs(value2 / value - 2.0 / 3.0) < 1.0e-12, "Perihelion advance does not scale as (2 + 2 gamma - beta)/3");
    value := GeneralRelativity.Perturbation.ppnLightDeflection(GeneralRelativity.Constants.GM_sun, 6.957e8) * arcsec;
    assert(abs(value / 1.75 - 1.0) < 5.0e-3, "Light deflection at the solar limb is not 1.75 arcsec");
    assert(abs(GeneralRelativity.Perturbation.ppnLightDeflection(GeneralRelativity.Constants.GM_sun, 6.957e8, 0.0) * arcsec
      / value - 0.5) < 1.0e-12, "gamma = 0 does not halve the deflection");
    // Shapiro delay near superior conjunction: (1 + gamma) GM/c^3 ln(4 r1 r2 / b^2).
    b := 6.957e8;
    r1 := 1.496e11;
    r2 := 2.279e11;
    chord := sqrt(r1^2 - b^2) + sqrt(r2^2 - b^2);
    value := GeneralRelativity.Perturbation.ppnShapiroDelay(GeneralRelativity.Constants.GM_sun, r1, r2, chord);
    value2 := 2.0 * GeneralRelativity.Constants.GM_sun / GeneralRelativity.Constants.c^3 * log(4.0 * r1 * r2 / b^2);
    assert(abs(value / value2 - 1.0) < 1.0e-4, "Shapiro delay does not reduce to (1 + gamma) GM/c^3 ln(4 r1 r2 / b^2)");
    assert(abs(2.0 * value / 247.2e-6 - 1.0) < 1.0e-2, "Earth-Mars Shapiro round trip at the solar limb is not about 247 microseconds");

    // ---------------------------------------------------------------- Gauge transformations
    // Pure gauge from a cubic xi_b: analytic ddh agrees with differenced h, and its linear Riemann vanishes.
    x := {0.4, -0.3, 0.7, 1.1};
    dddxi := gaugeThirdDerivative();
    for i in 1:4 loop
      for j in 1:4 loop
        hNew := GeneralRelativity.Perturbation.gaugeTransform(zeros(4, 4), dddxi[i, j, :, :]);
        for a in 1:4 loop
          for b in 1:4 loop
            ddh[i, j, a, b] := hNew[a, b];
          end for;
        end for;
      end for;
    end for;
    ddhFD := pureGaugeSecondDifference(x, 0.05);
    err := 0.0;
    for i in 1:4 loop
      for j in 1:4 loop
        err := max(err, Tests.Assertions.maxAbsMatrix(ddhFD[i, j, :, :] - ddh[i, j, :, :]));
      end for;
    end for;
    assert(err < 1.0e-9, "Differenced pure-gauge h does not match the analytic third derivatives of xi");
    R := GeneralRelativity.Perturbation.linearizedRiemann(ddh);
    R2 := GeneralRelativity.Perturbation.linearizedRiemann(ddhFD);
    err := 0.0;
    value := 0.0;
    for i in 1:4 loop
      for j in 1:4 loop
        err := max(err, Tests.Assertions.maxAbsMatrix(R[i, j, :, :]));
        err := max(err, Tests.Assertions.maxAbsMatrix(R2[i, j, :, :]) * 1.0e-4);
        value := max(value, Tests.Assertions.maxAbsMatrix(ddh[i, j, :, :]));
      end for;
    end for;
    assert(value > 0.1 and err < 1.0e-13, "Linearized Riemann tensor of a pure-gauge perturbation is not zero");
    // A non-gauge perturbation (static Newtonian potential) has curvature.
    hessian := [1.0, 0.2, -0.3; 0.2, -2.0, 0.4; -0.3, 0.4, 0.5];
    (R, Glin, laplacian) := GeneralRelativity.Perturbation.weakFieldCurvature(hessian, 1.0);
    assert(abs(R[1, 2, 1, 2] - hessian[1, 1]) < 1.0e-14, "Newtonian perturbation has no tidal curvature");
    // Killing vectors of Minkowski space (antisymmetric d_a xi_b) leave h unchanged.
    h := [0.2, 0.1, 0.0, 0.0; 0.1, 0.3, 0.05, 0.0; 0.0, 0.05, -0.1, 0.0; 0.0, 0.0, 0.0, 0.4];
    omega := [0.0, 0.3, -0.1, 0.2; -0.3, 0.0, 0.5, -0.4; 0.1, -0.5, 0.0, 0.7; -0.2, 0.4, -0.7, 0.0];
    assert(Tests.Assertions.maxAbsMatrix(GeneralRelativity.Perturbation.gaugeTransform(h, omega) - h) < 1.0e-15,
      "A Lorentz-transformation Killing vector changed h");
    assert(Tests.Assertions.maxAbsMatrix(GeneralRelativity.Perturbation.gaugeTransform(h, gaugeGradient(x))
      - (h - gaugeGradient(x) - transpose(gaugeGradient(x)))) < 1.0e-15, "gaugeTransform is not h - d_a xi_b - d_b xi_a");
    // TT plane wave h_+ = 0.01 cos(K.x), h_x = 0.02 cos(K.x) along z: curvature unchanged by the cubic gauge transformation.
    waveK := {-1.3, 0.0, 0.0, 1.3};
    phase := waveK * x;
    ddh := zeros(4, 4, 4, 4);
    for i in 1:4 loop
      for j in 1:4 loop
        ddh[i, j, 2, 2] := -0.01 * waveK[i] * waveK[j] * cos(phase);
        ddh[i, j, 3, 3] := 0.01 * waveK[i] * waveK[j] * cos(phase);
        ddh[i, j, 2, 3] := -0.02 * waveK[i] * waveK[j] * cos(phase);
        ddh[i, j, 3, 2] := ddh[i, j, 2, 3];
        hNew := GeneralRelativity.Perturbation.gaugeTransform(ddh[i, j, :, :], dddxi[i, j, :, :]);
        for a in 1:4 loop
          for b in 1:4 loop
            ddhNew[i, j, a, b] := hNew[a, b];
          end for;
        end for;
      end for;
    end for;
    R := GeneralRelativity.Perturbation.linearizedRiemann(ddh);
    R2 := GeneralRelativity.Perturbation.linearizedRiemann(ddhNew);
    err := 0.0;
    for i in 1:4 loop
      for j in 1:4 loop
        err := max(err, Tests.Assertions.maxAbsMatrix(R2[i, j, :, :] - R[i, j, :, :]));
      end for;
    end for;
    // R_0x0x = -1/2 d_t^2 h_xx for a TT wave.
    assert(abs(R[1, 2, 1, 2] - 0.5 * 0.01 * waveK[1]^2 * cos(phase)) < 1.0e-15, "TT wave R_0x0x is not -h_xx,tt / 2");
    assert(err < 1.0e-13, "Gauge transformation changed the linearized Riemann tensor of a TT wave");

    // ---------------------------------------------------------------- Electromagnetic stress-energy
    E := {0.3, -0.2, 0.5};
    B := {0.1, 0.4, -0.25};
    F := GeneralRelativity.StressEnergy.fieldStrength(E, B);
    T := GeneralRelativity.StressEnergy.electromagnetic(eta, F);
    assert(abs(GeneralRelativity.StressEnergy.trace(eta, T)) < 1.0e-15, "Electromagnetic T is not traceless");
    assert(abs(T[1, 1] - (E * E + B * B) / (8.0 * pi)) < 1.0e-15, "Gaussian energy density is not (E^2 + B^2) / (8 pi)");
    assert(Tests.Assertions.maxAbsVector(T[1, 2:4] + cross(E, B) / (4.0 * pi)) < 1.0e-15,
      "Gaussian T_0i is not minus the Poynting vector E x B / (4 pi)");
    // SI: T_00 = eps0 E^2 / 2 + B^2 / (2 mu0), T_0i = -S_i / c.
    cLight := GeneralRelativity.Constants.c;
    mu0 := 4.0e-7 * pi;
    E := {1.0e3, -2.0e2, 5.0e2};
    B := {1.0e-6, 3.0e-6, -2.0e-6};
    F := GeneralRelativity.StressEnergy.fieldStrength(E, B, cLight);
    T := GeneralRelativity.StressEnergy.electromagnetic(eta, F, mu0);
    value := 0.5 * (E * E) / (mu0 * cLight^2) + 0.5 * (B * B) / mu0;
    assert(abs(T[1, 1] / value - 1.0) < 1.0e-12, "SI electromagnetic energy density is not eps0 E^2/2 + B^2/(2 mu0)");
    assert(Tests.Assertions.maxAbsVector(T[1, 2:4] + cross(E, B) / (mu0 * cLight)) < 1.0e-12 * value,
      "SI electromagnetic T_0i is not -S/c");

    // ---------------------------------------------------------------- Polytrope equation of state
    value := GeneralRelativity.StressEnergy.polytropePressure(2.0, 3.0, 1.5);
    assert(abs(value - 3.0 * 2.0^1.5) < 1.0e-14, "polytropePressure is not K rho^Gamma");
    assert(abs(GeneralRelativity.StressEnergy.polytropeDensity(value, 3.0, 1.5) - 2.0) < 1.0e-14,
      "polytropeDensity does not invert polytropePressure");
    // Adiabatic first law de/drho = (e + p)/rho.
    dvalue := (GeneralRelativity.StressEnergy.polytropeEnergyDensity(2.0 + 1.0e-5, 3.0, 1.5, 2.0)
      - GeneralRelativity.StressEnergy.polytropeEnergyDensity(2.0 - 1.0e-5, 3.0, 1.5, 2.0)) / 2.0e-5;
    value2 := GeneralRelativity.StressEnergy.polytropeEnergyDensity(2.0, 3.0, 1.5, 2.0);
    assert(abs(dvalue - (value2 + value) / 2.0) < 1.0e-8, "Polytrope energy density violates de/drho = (e + p)/rho");
    assert(abs(value2 - (2.0 * 2.0^2 + value / 0.5)) < 1.0e-13, "Polytrope energy density is not rho c^2 + p/(Gamma - 1)");

    // ---------------------------------------------------------------- Harmonics: quadrature and Legendre functions
    (nodes, weights) := GeneralRelativity.Harmonics.gaussLegendre(6);
    assert(abs(sum(weights) - 2.0) < 1.0e-14 and abs(weights * nodes .^ 10 - 2.0 / 11.0) < 1.0e-14
      and abs(weights * nodes .^ 11) < 1.0e-14, "6-point Gauss-Legendre does not integrate degree 11 exactly");
    assert(nodes[1] < nodes[2] and abs(nodes[1] + nodes[6]) < 1.0e-15, "Gauss-Legendre nodes are not ascending and symmetric");
    value := 0.3;
    assert(abs(GeneralRelativity.Harmonics.associatedLegendre(2, 1, value) + 3.0 * value * sqrt(1.0 - value^2)) < 1.0e-15
      and abs(GeneralRelativity.Harmonics.associatedLegendre(2, 2, value) - 3.0 * (1.0 - value^2)) < 1.0e-15
      and abs(GeneralRelativity.Harmonics.associatedLegendre(3, 0, value) - 0.5 * (5.0 * value^3 - 3.0 * value)) < 1.0e-15
      and abs(GeneralRelativity.Harmonics.associatedLegendre(1, 1, value) + sqrt(1.0 - value^2)) < 1.0e-15,
      "Associated Legendre functions miss their closed forms (Condon-Shortley phase)");
    assert(abs(GeneralRelativity.Harmonics.associatedLegendre(2, -1, value)
      + GeneralRelativity.Harmonics.associatedLegendre(2, 1, value) / 6.0) < 1.0e-15,
      "P_l^-m is not (-1)^m (l - m)!/(l + m)! P_l^m");
    assert(abs(GeneralRelativity.Harmonics.normalizedLegendre(7, -3, value)
      - sqrt(15.0 / (4.0 * pi) * exp(GeneralRelativity.Harmonics.logFactorial(10) - GeneralRelativity.Harmonics.logFactorial(4)))
        * GeneralRelativity.Harmonics.associatedLegendre(7, -3, value)) < 1.0e-13,
      "normalizedLegendre is not the normalized associatedLegendre for negative m");
    assert(abs(GeneralRelativity.Harmonics.logFactorial(10) - log(3628800.0)) < 1.0e-13, "logFactorial(10) is wrong");
    // Y_10 and Y_11 closed forms; conjugation symmetry.
    theta := 0.7;
    phi := 1.9;
    (re, im) := GeneralRelativity.Harmonics.sphericalHarmonic(1, 0, theta, phi);
    assert(abs(re - sqrt(3.0 / (4.0 * pi)) * cos(theta)) < 1.0e-15 and abs(im) < 1.0e-15, "Y_10 is not sqrt(3/(4 pi)) cos(theta)");
    (re, im) := GeneralRelativity.Harmonics.sphericalHarmonic(1, 1, theta, phi);
    assert(abs(re + sqrt(3.0 / (8.0 * pi)) * sin(theta) * cos(phi)) < 1.0e-15
      and abs(im + sqrt(3.0 / (8.0 * pi)) * sin(theta) * sin(phi)) < 1.0e-15, "Y_11 is not -sqrt(3/(8 pi)) sin(theta) e^(i phi)");
    (re, im) := GeneralRelativity.Harmonics.sphericalHarmonic(4, 3, theta, phi);
    (re2, im2) := GeneralRelativity.Harmonics.sphericalHarmonic(4, -3, theta, phi);
    assert(abs(re2 + re) < 1.0e-15 and abs(im2 - im) < 1.0e-15, "Y_l,-m is not (-1)^m conj(Y_lm)");
    assert(abs(GeneralRelativity.Harmonics.realSphericalHarmonic(1, 1, theta, phi) - sqrt(3.0 / (4.0 * pi)) * sin(theta) * cos(phi)) < 1.0e-15
      and abs(GeneralRelativity.Harmonics.realSphericalHarmonic(1, -1, theta, phi) - sqrt(3.0 / (4.0 * pi)) * sin(theta) * sin(phi)) < 1.0e-15,
      "Real harmonics Y_1,1 and Y_1,-1 are not sqrt(3/(4 pi)) x/r and y/r");
    // Addition theorem sum_m |Y_lm|^2 = (2l + 1)/(4 pi), up to l = 20.
    for l in {2, 4, 10, 20} loop
      sumSq := 0.0;
      for m in -l:l loop
        (re, im) := GeneralRelativity.Harmonics.sphericalHarmonic(l, m, theta, phi);
        sumSq := sumSq + re^2 + im^2;
      end for;
      assert(abs(sumSq / ((2 * l + 1) / (4.0 * pi)) - 1.0) < 1.0e-12, "Addition theorem fails at l = " + String(l));
    end for;
    // Orthonormality of complex and real harmonics for l <= 4: 6 Gauss-Legendre nodes in cos(theta) x 12 azimuths.
    idx := 0;
    for l in 0:4 loop
      for m in -l:l loop
        idx := idx + 1;
        lmTable[idx, :] := {l, m};
      end for;
    end for;
    for i in 1:6 loop
      for j in 1:12 loop
        theta := acos(nodes[i]);
        phi := 2.0 * pi * (j - 1) / 12.0;
        w[(i - 1) * 12 + j] := weights[i] * 2.0 * pi / 12.0;
        for n in 1:25 loop
          (re, im) := GeneralRelativity.Harmonics.sphericalHarmonic(lmTable[n, 1], lmTable[n, 2], theta, phi);
          Yre[n, (i - 1) * 12 + j] := re;
          Yim[n, (i - 1) * 12 + j] := im;
          Yreal[n, (i - 1) * 12 + j] := GeneralRelativity.Harmonics.realSphericalHarmonic(lmTable[n, 1], lmTable[n, 2], theta, phi);
        end for;
      end for;
    end for;
    err := 0.0;
    for n in 1:25 loop
      for n2 in 1:25 loop
        gram := sum(w[q] * (Yre[n, q] * Yre[n2, q] + Yim[n, q] * Yim[n2, q]) for q in 1:72);
        quadIm := sum(w[q] * (Yim[n, q] * Yre[n2, q] - Yre[n, q] * Yim[n2, q]) for q in 1:72);
        gramReal := sum(w[q] * Yreal[n, q] * Yreal[n2, q] for q in 1:72);
        err := max(err, max(abs(gram - (if n == n2 then 1.0 else 0.0)), max(abs(quadIm), abs(gramReal - (if n == n2 then 1.0 else 0.0)))));
      end for;
    end for;
    assert(err < 1.0e-13, "Spherical harmonics with l <= 4 are not orthonormal");

    // ---------------------------------------------------------------- Clebsch-Gordan and 3j symbols
    assert(abs(GeneralRelativity.Harmonics.clebschGordan(0.5, 0.5, 0.5, -0.5, 0.0, 0.0) - 1.0 / sqrt(2.0)) < 1.0e-15
      and abs(GeneralRelativity.Harmonics.clebschGordan(0.5, -0.5, 0.5, 0.5, 0.0, 0.0) + 1.0 / sqrt(2.0)) < 1.0e-15,
      "Singlet Clebsch-Gordan coefficients are not +-1/sqrt(2)");
    assert(abs(GeneralRelativity.Harmonics.clebschGordan(1, 0, 1, 0, 2, 0) - sqrt(2.0 / 3.0)) < 1.0e-15
      and abs(GeneralRelativity.Harmonics.clebschGordan(1, 1, 1, -1, 0, 0) - 1.0 / sqrt(3.0)) < 1.0e-15
      and abs(GeneralRelativity.Harmonics.clebschGordan(1, 0, 1, 0, 1, 0)) < 1.0e-15
      and abs(GeneralRelativity.Harmonics.clebschGordan(1, 1, 1, 0, 1, 0)) < 1.0e-15,
      "Clebsch-Gordan coefficients for 1 x 1 are wrong");
    assert(abs(GeneralRelativity.Harmonics.clebschGordan(7, 7, 6.5, 6.5, 13.5, 13.5) - 1.0) < 1.0e-13,
      "Stretched Clebsch-Gordan coefficient is not 1");
    // Orthogonality sum_{m1,m2} <j1 m1 j2 m2|J M><j1 m1 j2 m2|J' M'> = delta for j1 = 1, j2 = 3/2.
    err := 0.0;
    for J in {0.5, 1.5, 2.5} loop
      for J2 in {0.5, 1.5, 2.5} loop
        for M in -2.5:2.5 loop
          value := 0.0;
          for m1 in -1:1 loop
            for m2 in {-1.5, -0.5, 0.5, 1.5} loop
              value := value + GeneralRelativity.Harmonics.clebschGordan(1, m1, 1.5, m2, J, M)
                * GeneralRelativity.Harmonics.clebschGordan(1, m1, 1.5, m2, J2, M);
            end for;
          end for;
          err := max(err, abs(value - (if J == J2 and abs(M) <= J then 1.0 else 0.0)));
        end for;
      end for;
    end for;
    assert(err < 1.0e-14, "Clebsch-Gordan coefficients for 1 x 3/2 are not orthonormal");
    // Orthogonality at high angular momentum, j1 = 10, j2 = 12, M = 1, all J.
    err := 0.0;
    for J in 2:22 loop
      for J2 in 2:22 loop
        value := sum(GeneralRelativity.Harmonics.clebschGordan(10, m1, 12, 1 - m1, J, 1)
          * GeneralRelativity.Harmonics.clebschGordan(10, m1, 12, 1 - m1, J2, 1) for m1 in -10:10);
        err := max(err, abs(value - (if J == J2 then 1.0 else 0.0)));
      end for;
    end for;
    assert(err < 1.0e-10, "Clebsch-Gordan coefficients for 10 x 12 are not orthonormal");
    // 3j symbols: closed forms and symmetries.
    assert(abs(GeneralRelativity.Harmonics.wigner3j(1, 1, 2, 0, 0, 0) - sqrt(2.0 / 15.0)) < 1.0e-15
      and abs(GeneralRelativity.Harmonics.wigner3j(2, 2, 2, 0, 0, 0) + sqrt(2.0 / 35.0)) < 1.0e-15
      and abs(GeneralRelativity.Harmonics.wigner3j(1, 1, 1, 0, 0, 0)) < 1.0e-15,
      "3j symbols (1 1 2; 0 0 0), (2 2 2; 0 0 0) or (1 1 1; 0 0 0) are wrong");
    w1 := GeneralRelativity.Harmonics.wigner3j(2, 3, 4, 1, -2, 1);
    assert(abs(w1) > 1.0e-2, "Test 3j symbol unexpectedly vanishes");
    assert(abs(GeneralRelativity.Harmonics.wigner3j(3, 4, 2, -2, 1, 1) - w1) < 1.0e-15
      and abs(GeneralRelativity.Harmonics.wigner3j(4, 2, 3, 1, 1, -2) - w1) < 1.0e-15,
      "3j symbol is not invariant under cyclic permutations");
    assert(abs(GeneralRelativity.Harmonics.wigner3j(3, 2, 4, -2, 1, 1) + w1) < 1.0e-15
      and abs(GeneralRelativity.Harmonics.wigner3j(2, 3, 4, -1, 2, -1) + w1) < 1.0e-15,
      "3j symbol does not change by (-1)^(j1 + j2 + j3) = -1 under a swap or m -> -m");
    w2 := GeneralRelativity.Harmonics.wigner3j(1.5, 2, 2.5, 0.5, -1, 0.5);
    assert(abs(w2) > 1.0e-2 and abs(GeneralRelativity.Harmonics.wigner3j(2, 1.5, 2.5, -1, 0.5, 0.5) - w2) < 1.0e-15,
      "Half-integer 3j symbol has the wrong swap phase (j1 + j2 + j3 = 6 is even)");

    // ---------------------------------------------------------------- Gaunt coefficients against quadrature
    (nodes10, weights10) := GeneralRelativity.Harmonics.gaussLegendre(10);
    gauntCases := [1, 1, 1, -1, 2, 0; 2, 1, 3, -2, 3, 1; 2, 2, 2, -1, 4, -1; 1, 0, 1, 0, 1, 0];
    for c in 1:4 loop
      quadRe := 0.0;
      quadIm := 0.0;
      for i in 1:10 loop
        for j in 1:16 loop
          theta := acos(nodes10[i]);
          phi := 2.0 * pi * (j - 1) / 16.0;
          (re, im) := GeneralRelativity.Harmonics.sphericalHarmonic(gauntCases[c, 1], gauntCases[c, 2], theta, phi);
          (re2, im2) := GeneralRelativity.Harmonics.sphericalHarmonic(gauntCases[c, 3], gauntCases[c, 4], theta, phi);
          (re3, im3) := GeneralRelativity.Harmonics.sphericalHarmonic(gauntCases[c, 5], gauntCases[c, 6], theta, phi);
          // Real and imaginary parts of the triple product.
          quadRe := quadRe + weights10[i] * 2.0 * pi / 16.0 * ((re * re2 - im * im2) * re3 - (re * im2 + im * re2) * im3);
          quadIm := quadIm + weights10[i] * 2.0 * pi / 16.0 * ((re * re2 - im * im2) * im3 + (re * im2 + im * re2) * re3);
        end for;
      end for;
      value := GeneralRelativity.Harmonics.gaunt(gauntCases[c, 1], gauntCases[c, 2], gauntCases[c, 3], gauntCases[c, 4],
        gauntCases[c, 5], gauntCases[c, 6]);
      assert(abs(value - quadRe) < 1.0e-14 and abs(quadIm) < 1.0e-14, "Gaunt coefficient disagrees with quadrature in case " + String(c));
      assert(c == 4 or abs(value) > 1.0e-2, "Nonzero Gaunt test case vanished: " + String(c));
    end for;
    assert(abs(GeneralRelativity.Harmonics.gaunt(1, 0, 1, 0, 1, 0)) < 1.0e-15
      and abs(GeneralRelativity.Harmonics.gaunt(1, 1, 1, 0, 2, 0)) < 1.0e-15
      and abs(GeneralRelativity.Harmonics.gaunt(1, 0, 1, 0, 4, 0)) < 1.0e-15,
      "Gaunt coefficient violates parity, m-sum or triangle selection rules");
    assert(abs(GeneralRelativity.Harmonics.gaunt(0, 0, 3, 2, 3, -2) - 1.0 / sqrt(4.0 * pi)) < 1.0e-15,
      "Gaunt coefficient with Y_00 is not (-1)^m / sqrt(4 pi)");

    passed := true;
  end run;

  parameter Boolean passed = run();
equation
  assert(passed, "GeneralRelativity metric, source, PPN, gauge and harmonics assertions did not complete");
  annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 1.0e-8, Interval = 0.1));
end GeneralRelativityMetricTests;
