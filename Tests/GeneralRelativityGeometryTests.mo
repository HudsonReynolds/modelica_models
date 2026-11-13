within Tests;
model GeneralRelativityGeometryTests "3+1 foliation and constraints, Killing vectors and Lie derivatives, exterior calculus and Cartan structure equations"
  // ---------------------------------------------------------------------------
  // Test-local metrics.
  function deSitterFlat "Flat-slicing de Sitter metric diag(-1, a^2, a^2, a^2), a = exp(H x0); p = {H}"
    extends GeneralRelativity.Interfaces.MetricFunction;
  protected
    Real a2 = exp(2.0 * p[1] * x[1]);
  algorithm
    g := diagonal({-1.0, a2, a2, a2});
  end deSitterFlat;

  function painleveGullstrand "Schwarzschild in ingoing Painleve-Gullstrand coordinates {cT, r, theta, phi}; p = {rs}"
    extends GeneralRelativity.Interfaces.MetricFunction;
  protected
    Real r = x[2];
  algorithm
    g := zeros(4, 4);
    g[1, 1] := -(1.0 - p[1] / r);
    g[1, 2] := sqrt(p[1] / r);
    g[2, 1] := g[1, 2];
    g[2, 2] := 1.0;
    g[3, 3] := r^2;
    g[4, 4] := r^2 * sin(x[3])^2;
  end painleveGullstrand;

  function wobblyMetric "Time-dependent, shifted, non-vacuum Lorentzian metric in {ct, x, y, z} with no symmetry"
    extends GeneralRelativity.Interfaces.MetricFunction;
  protected
    Real t = x[1];
    Real X = x[2];
    Real Y = x[3];
    Real Z = x[4];
  algorithm
    g[1, 1] := -(1.0 + 0.3 * sin(t + X)^2 + 0.1 * Y^2);
    g[1, 2] := 0.1 * sin(Y + t);
    g[1, 3] := 0.05 * X * Z;
    g[1, 4] := 0.1 * cos(X) * t;
    g[2, 2] := 1.0 + 0.1 * X^2 + 0.05 * t;
    g[2, 3] := 0.1 * sin(Z * t);
    g[2, 4] := 0.02 * X * Y;
    g[3, 3] := 1.0 + 0.2 * cos(Y) * exp(0.1 * t);
    g[3, 4] := 0.03 * t * Z;
    g[4, 4] := 1.5 + 0.1 * sin(X + Z);
    for i in 2:4 loop
      for j in 1:i - 1 loop
        g[i, j] := g[j, i];
      end for;
    end for;
  end wobblyMetric;

  function reissnerNordstrom "Reissner-Nordstrom metric in {ct, r, theta, phi}; p = {rs, rQ}, f = 1 - rs/r + rQ^2/r^2"
    extends GeneralRelativity.Interfaces.MetricFunction;
  protected
    Real r = x[2];
    Real f = 1.0 - p[1] / r + p[2]^2 / r^2;
  algorithm
    g := diagonal({-f, 1.0 / f, r^2, r^2 * sin(x[3])^2});
  end reissnerNordstrom;

  // ---------------------------------------------------------------------------
  // Test-local forms and tetrads.
  function radialFieldStrength "F = (Q / r^n) dx0 ^ dr; p = {Q, n}; n = 2 is the Coulomb field"
    extends GeneralRelativity.ExteriorCalculus.TwoFormField;
  algorithm
    F := zeros(4, 4);
    F[1, 2] := p[1] / x[2]^p[2];
    F[2, 1] := -F[1, 2];
  end radialFieldStrength;

  function testScalar "Smooth scalar field"
    extends GeneralRelativity.ExteriorCalculus.ScalarField;
  algorithm
    f := sin(x[1] * x[2]) + x[3]^2 * cos(x[4]) + exp(0.2 * x[2]) * x[4];
  end testScalar;

  function gradientOfTestScalar "d f of testScalar by central differences"
    extends GeneralRelativity.ExteriorCalculus.OneFormField;
  algorithm
    w := GeneralRelativity.ExteriorCalculus.exteriorDerivative0(testScalar, x, p, {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3});
  end gradientOfTestScalar;

  function testOneForm "Smooth 1-form"
    extends GeneralRelativity.ExteriorCalculus.OneFormField;
  algorithm
    w := {x[3] * sin(x[2]), x[1]^2 + cos(x[4]), x[2] * x[4], exp(0.1 * x[1]) * x[3]};
  end testOneForm;

  function positionOneForm "Lowered Minkowski position x_a = eta_ab x^b"
    extends GeneralRelativity.ExteriorCalculus.OneFormField;
  algorithm
    w := {-x[1], x[2], x[3], x[4]};
  end positionOneForm;

  function dOfTestOneForm "d A of testOneForm by central differences"
    extends GeneralRelativity.ExteriorCalculus.TwoFormField;
  algorithm
    F := GeneralRelativity.ExteriorCalculus.exteriorDerivative1(testOneForm, x, p, {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3});
  end dOfTestOneForm;

  function testTwoForm "Smooth non-simple 2-form"
    extends GeneralRelativity.ExteriorCalculus.TwoFormField;
  protected
    Real a[4] = {x[2], 0.5, sin(x[3]), x[1] * x[4]};
    Real b[4] = {cos(x[4]), x[1], 0.3 * x[2], 1.0};
    Real c[4] = {0.2, x[3]^2, 1.0, -x[2]};
    Real d[4] = {x[4], 0.1, x[1], 0.7};
  algorithm
    F := GeneralRelativity.ExteriorCalculus.wedge11(a, b) + GeneralRelativity.ExteriorCalculus.wedge11(c, d);
  end testTwoForm;

  function dOfTestTwoForm "d F of testTwoForm by central differences"
    extends GeneralRelativity.ExteriorCalculus.ThreeFormField;
  algorithm
    H := GeneralRelativity.ExteriorCalculus.exteriorDerivative2(testTwoForm, x, p, {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3});
  end dOfTestTwoForm;

  function dualOfTestOneForm "Schwarzschild Hodge dual of testOneForm; p = {rs}"
    extends GeneralRelativity.ExteriorCalculus.ThreeFormField;
  protected
    Real w[4];
    Real g[4, 4];
  algorithm
    w := testOneForm(x, p);
    g := GeneralRelativity.Metrics.schwarzschild(x, p);
    H := GeneralRelativity.ExteriorCalculus.hodge1(w, g);
  end dualOfTestOneForm;

  function dualOfTestTwoForm "Schwarzschild Hodge dual of testTwoForm; p = {rs}"
    extends GeneralRelativity.ExteriorCalculus.TwoFormField;
  protected
    Real G2[4, 4];
    Real g[4, 4];
  algorithm
    G2 := testTwoForm(x, p);
    g := GeneralRelativity.Metrics.schwarzschild(x, p);
    F := GeneralRelativity.ExteriorCalculus.hodge2(G2, g);
  end dualOfTestTwoForm;

  function rotatedTetrad "Static Schwarzschild tetrad after a position-dependent boost (e0, e1) and rotation (e2, e3); p = {rs}"
    extends GeneralRelativity.ExteriorCalculus.TetradField;
  protected
    Real s[4, 4];
    Real rapidity = 0.3 * sin(x[3]) + 0.05 * x[2] + 0.1 * x[1];
    Real angle = 0.5 * x[4] + 0.2 * cos(x[3]) + 0.03 * x[2];
  algorithm
    s := GeneralRelativity.ExteriorCalculus.schwarzschildTetrad(x, p);
    e[1, :] := cosh(rapidity) * s[1, :] + sinh(rapidity) * s[2, :];
    e[2, :] := sinh(rapidity) * s[1, :] + cosh(rapidity) * s[2, :];
    e[3, :] := cos(angle) * s[3, :] + sin(angle) * s[4, :];
    e[4, :] := -sin(angle) * s[3, :] + cos(angle) * s[4, :];
  end rotatedTetrad;

  // ---------------------------------------------------------------------------
  function maxAbs3 "Largest absolute element of a rank-3 array"
    input Real A[:, :, :];
    output Real m;
  algorithm
    m := max(abs(A[i, j, k]) for i in 1:size(A, 1), j in 1:size(A, 2), k in 1:size(A, 3));
  end maxAbs3;

  function maxAbs4 "Largest absolute element of a rank-4 array"
    input Real A[:, :, :, :];
    output Real m;
  algorithm
    m := max(abs(A[i, j, k, l]) for i in 1:size(A, 1), j in 1:size(A, 2), k in 1:size(A, 3), l in 1:size(A, 4));
  end maxAbs4;

  // ---------------------------------------------------------------------------
  function runFoliation "3+1 split, extrinsic curvature and Gauss-Codazzi constraints"
    output Boolean passed;
  protected
    constant Real rs = 2.0;
    constant Real step[4] = {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3};
    GeneralRelativity.Foliation.ADM s;
    GeneralRelativity.Tensors.Curvature k;
    Real x[4];
    Real g[4, 4];
    Real dg[4, 4, 4];
    Real ddg[4, 4, 4, 4];
    Real dK[4, 3, 3];
    Real dKNumeric[3, 3];
    Real KPlus[3, 3];
    Real KMinus[3, 3];
    Real xPlus[4];
    Real xMinus[4];
    Real nUp[4];
    Real nDown[4];
    Real alpha;
    Real betaUp[3];
    Real betaDown[3];
    Real gamma[3, 3];
    Real gammaInv[3, 3];
    Boolean accepted;
    Real KUp[3, 3];
    Real scale;
    Real KExact[3, 3];
    Real H;
    Real r;
    Real a2;
  algorithm
    // Algebra of the split on a generic metric.
    x := {0.3, 0.4, -0.2, 0.5};
    g := wobblyMetric(x, {0.0});
    (alpha, betaUp, betaDown, gamma, gammaInv, accepted) := GeneralRelativity.Foliation.admDecomposition(g);
    assert(accepted, "ADM split of the test metric was rejected");
    assert(abs(-alpha^2 + betaUp * betaDown - g[1, 1]) < 1.0e-12
      and Tests.Assertions.maxAbsVector(betaDown - g[1, 2:4]) < 1.0e-15
      and Tests.Assertions.maxAbsMatrix(gamma * gammaInv - identity(3)) < 1.0e-12,
      "ADM split does not reconstruct the metric");
    (nUp, nDown) := GeneralRelativity.Foliation.normalVector(g);
    assert(abs(nUp * g * nUp + 1.0) < 1.0e-12 and Tests.Assertions.maxAbsVector(g * nUp - nDown) < 1.0e-12,
      "Normal vector is not unit or not the dual of n_a = (-alpha, 0, 0, 0)");
    assert(nUp[1] > 0.0, "Normal vector is not future-pointing");

    // Generic time-dependent, shifted, non-vacuum metric: K_ij = -alpha Gamma^0_ij and both constraints hold.
    s := GeneralRelativity.Foliation.admAt(wobblyMetric, x, {0.0}, step);
    k := GeneralRelativity.Metrics.curvatureAt(wobblyMetric, x, {0.0}, step);
    assert(s.valid > 0.5, "ADM record of the test metric is not valid");
    assert(Tests.Assertions.maxAbsMatrix(s.K + s.alpha * k.Gamma[1, 2:4, 2:4]) < 1.0e-10,
      "Extrinsic curvature disagrees with -alpha Gamma^0_ij");
    assert(Tests.Assertions.maxAbsMatrix(s.K) > 1.0e-2, "Test metric has negligible extrinsic curvature; the test is too weak");
    assert(abs(s.normalEinstein) > 1.0e-2 and Tests.Assertions.maxAbsVector(s.momentumEinstein) > 1.0e-2,
      "Test metric is nearly vacuum; the constraint test is too weak");
    KUp := s.KMixed * s.gammaInv;
    scale := abs(s.ricciScalar3) + s.trK^2 + abs(sum(s.K[i, j] * KUp[i, j] for i in 1:3, j in 1:3)) + 2.0 * abs(s.normalEinstein);
    assert(abs(s.hamiltonian) < 1.0e-9 * scale, "Gauss (Hamiltonian) constraint residual is not zero: " + String(s.hamiltonian));
    assert(Tests.Assertions.maxAbsVector(s.momentum) < 1.0e-9 * (1.0 + Tests.Assertions.maxAbsVector(s.momentumEinstein)),
      "Codazzi (momentum) constraint residual is not zero: " + String(Tests.Assertions.maxAbsVector(s.momentum)));
    // The chain-rule gradient of K agrees with nested differencing (outer step 2e-3, inner 1e-3).
    g := wobblyMetric(x, {0.0});
    dg := GeneralRelativity.Metrics.derivative(wobblyMetric, x, {0.0}, step);
    ddg := GeneralRelativity.Metrics.secondDerivative(wobblyMetric, x, {0.0}, step);
    dK := GeneralRelativity.Foliation.extrinsicCurvatureGradient(g, dg, ddg);
    for kk in 1:4 loop
      xPlus := x;
      xMinus := x;
      xPlus[kk] := x[kk] + 2.0e-3;
      xMinus[kk] := x[kk] - 2.0e-3;
      g := wobblyMetric(xPlus, {0.0});
      dg := GeneralRelativity.Metrics.derivative(wobblyMetric, xPlus, {0.0}, step);
      KPlus := GeneralRelativity.Foliation.extrinsicCurvature(g, dg);
      g := wobblyMetric(xMinus, {0.0});
      dg := GeneralRelativity.Metrics.derivative(wobblyMetric, xMinus, {0.0}, step);
      KMinus := GeneralRelativity.Foliation.extrinsicCurvature(g, dg);
      dKNumeric := (KPlus - KMinus) / 4.0e-3;
      assert(Tests.Assertions.maxAbsMatrix(dK[kk, :, :] - dKNumeric) < 1.0e-5,
        "extrinsicCurvatureGradient disagrees with differenced K along coordinate " + String(kk));
    end for;

    // Flat de Sitter slicing: K_ij = -H a^2 delta_ij, K = -3H < 0 (expanding), flat slices, G_nn = 3 H^2.
    H := 0.3;
    x := {0.5, 1.0, 2.0, -1.0};
    a2 := exp(2.0 * H * x[1]);
    s := GeneralRelativity.Foliation.admAt(deSitterFlat, x, {H}, step);
    assert(abs(s.alpha - 1.0) < 1.0e-12 and Tests.Assertions.maxAbsVector(s.betaUp) < 1.0e-12, "de Sitter lapse or shift is wrong");
    assert(Tests.Assertions.maxAbsMatrix(s.K + H * a2 * identity(3)) < 1.0e-6, "de Sitter extrinsic curvature is not -H a^2 delta");
    assert(abs(s.trK + 3.0 * H) < 1.0e-6, "de Sitter mean curvature is not -3H (sign convention)");
    assert(abs(s.ricciScalar3) < 1.0e-6, "Flat de Sitter slices are curved");
    assert(abs(s.normalEinstein - 3.0 * H^2) < 1.0e-5, "de Sitter energy density G_nn is not 3 H^2");
    assert(abs(s.hamiltonian) < 1.0e-9 and Tests.Assertions.maxAbsVector(s.momentum) < 1.0e-9, "de Sitter constraint residuals are not zero");

    // Painleve-Gullstrand: unit lapse, flat slices, K_rr = -sqrt(rs)/(2 r^1.5), K_thth = sqrt(rs r).
    x := {0.2, 5.0 * rs, 1.1, 0.4};
    r := x[2];
    s := GeneralRelativity.Foliation.admAt(painleveGullstrand, x, {rs}, step);
    assert(abs(s.alpha - 1.0) < 1.0e-12 and abs(s.betaUp[1] - sqrt(rs / r)) < 1.0e-12, "Painleve-Gullstrand lapse or shift is wrong");
    KExact := diagonal({-0.5 * sqrt(rs) * r^(-1.5), sqrt(rs * r), sqrt(rs * r) * sin(x[3])^2});
    assert(Tests.Assertions.maxAbsMatrix(s.K - KExact) < 1.0e-6, "Painleve-Gullstrand extrinsic curvature is wrong");
    assert(abs(s.trK - 1.5 * sqrt(rs) * r^(-1.5)) < 1.0e-7, "Painleve-Gullstrand mean curvature is not +3/2 sqrt(rs) r^-1.5");
    assert(abs(s.ricciScalar3) < 1.0e-6 and abs(s.normalEinstein) < 1.0e-6, "Painleve-Gullstrand slices are not flat vacuum slices");
    assert(abs(s.hamiltonian) < 1.0e-9 and Tests.Assertions.maxAbsVector(s.momentum) < 1.0e-9, "Painleve-Gullstrand constraint residuals are not zero");

    // Schwarzschild time: a time-symmetric slice, K = 0 and (3)R = 0 although (3)R_ij is not zero.
    x := {0.0, 5.0 * rs, 1.1, 0.4};
    s := GeneralRelativity.Foliation.admAt(GeneralRelativity.Metrics.schwarzschild, x, {rs}, step);
    assert(abs(s.alpha - sqrt(1.0 - rs / x[2])) < 1.0e-12, "Schwarzschild lapse is not sqrt(1 - rs/r)");
    assert(Tests.Assertions.maxAbsMatrix(s.K) < 1.0e-12, "Schwarzschild static slice has extrinsic curvature");
    assert(abs(s.ricciScalar3) < 1.0e-6, "Schwarzschild static slice has nonzero scalar curvature");
    assert(abs(s.ricci3[1, 1] + rs / (x[2]^3 * (1.0 - rs / x[2]))) < 1.0e-6,
      "Schwarzschild slice Ricci (3)R_rr is not -rs/(r^3 f)");
    passed := true;
  end runFoliation;

  // ---------------------------------------------------------------------------
  function runSymmetry "Lie derivatives, Killing and conformal Killing vectors"
    output Boolean passed;
  protected
    constant Real rs = 2.0;
    constant Real step[4] = {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3};
    constant Real eta[4, 4] = diagonal({-1.0, 1.0, 1.0, 1.0});
    Real xs[2, 4] = {{0.3, 5.0 * rs, 1.1, 0.4}, {-1.0, 3.5 * rs, 0.6, 2.0}};
    Real x[4];
    Real L[4, 4];
    Real scale;
    Real f;
  algorithm
    for i in 1:2 loop
      x := xs[i, :];
      // Schwarzschild: time translation and the full rotation group.
      assert(GeneralRelativity.Symmetry.isKilling(GeneralRelativity.Metrics.schwarzschild, GeneralRelativity.Symmetry.timeTranslation, x, {rs}, {0.0}, step),
        "d_t is not Killing in Schwarzschild");
      assert(GeneralRelativity.Symmetry.isKilling(GeneralRelativity.Metrics.schwarzschild, GeneralRelativity.Symmetry.axialRotation, x, {rs}, {0.0}, step),
        "d_phi is not Killing in Schwarzschild");
      assert(GeneralRelativity.Symmetry.isKilling(GeneralRelativity.Metrics.schwarzschild, GeneralRelativity.Symmetry.rotationX, x, {rs}, {0.0}, step),
        "rotationX is not Killing in Schwarzschild");
      assert(GeneralRelativity.Symmetry.isKilling(GeneralRelativity.Metrics.schwarzschild, GeneralRelativity.Symmetry.rotationY, x, {rs}, {0.0}, step),
        "rotationY is not Killing in Schwarzschild");
      assert(GeneralRelativity.Symmetry.killingResidual(GeneralRelativity.Metrics.schwarzschild, GeneralRelativity.Symmetry.coordinateVector, x, {rs}, {2.0}, step) > 0.5,
        "d_r passed as a Killing vector of Schwarzschild");
      // Kerr keeps d_t and d_phi but loses the other rotations.
      assert(GeneralRelativity.Symmetry.isKilling(GeneralRelativity.Metrics.kerr, GeneralRelativity.Symmetry.timeTranslation, x, {rs, 0.4 * rs}, {0.0}, step),
        "d_t is not Killing in Kerr");
      assert(GeneralRelativity.Symmetry.isKilling(GeneralRelativity.Metrics.kerr, GeneralRelativity.Symmetry.axialRotation, x, {rs, 0.4 * rs}, {0.0}, step),
        "d_phi is not Killing in Kerr");
      assert(GeneralRelativity.Symmetry.killingResidual(GeneralRelativity.Metrics.kerr, GeneralRelativity.Symmetry.rotationX, x, {rs, 0.4 * rs}, {0.0}, step) > 1.0e-3,
        "rotationX passed as a Killing vector of Kerr");
      assert(GeneralRelativity.Symmetry.isKilling(GeneralRelativity.Metrics.kerr, GeneralRelativity.Symmetry.rotationX, x, {rs, 0.0}, {0.0}, step),
        "rotationX is not Killing in Kerr with a = 0");
    end for;

    // Closed forms in Schwarzschild: L_{d_r} g has (rr) = d_r g_rr; div d_r = 2/r.
    x := xs[1, :];
    f := 1.0 - rs / x[2];
    (L, scale) := GeneralRelativity.Symmetry.lieDerivativeMetric(GeneralRelativity.Metrics.schwarzschild, GeneralRelativity.Symmetry.coordinateVector, x, {rs}, {2.0}, step);
    assert(abs(L[2, 2] + rs / (x[2]^2 * f^2)) < 1.0e-6 and abs(L[3, 3] - 2.0 * x[2]) < 1.0e-6 and abs(L[1, 1] + rs / x[2]^2) < 1.0e-6,
      "Lie derivative of the Schwarzschild metric along d_r is wrong");
    assert(abs(GeneralRelativity.Symmetry.divergence(GeneralRelativity.Metrics.schwarzschild, GeneralRelativity.Symmetry.coordinateVector, x, {rs}, {2.0}, step) - 2.0 / x[2]) < 1.0e-7,
      "Divergence of d_r in Schwarzschild is not 2/r");
    assert(GeneralRelativity.Symmetry.conformalKillingResidual(GeneralRelativity.Metrics.schwarzschild, GeneralRelativity.Symmetry.coordinateVector, x, {rs}, {2.0}, step) > 0.1,
      "d_r passed as a conformal Killing vector of Schwarzschild");

    // Minkowski: the Poincare generators are Killing.
    x := {0.7, 1.3, -0.4, 2.1};
    assert(GeneralRelativity.Symmetry.isKilling(GeneralRelativity.Metrics.minkowski, GeneralRelativity.Symmetry.cartesianTranslation, x, {0.0}, {1.0, 2.0, 3.0, 4.0}, step),
      "Translations are not Killing in Minkowski");
    assert(GeneralRelativity.Symmetry.isKilling(GeneralRelativity.Metrics.minkowski, GeneralRelativity.Symmetry.cartesianRotation, x, {0.0}, {0.3, -0.5, 0.8}, step),
      "Rotations are not Killing in Minkowski");
    assert(GeneralRelativity.Symmetry.isKilling(GeneralRelativity.Metrics.minkowski, GeneralRelativity.Symmetry.cartesianBoost, x, {0.0}, {0.6, 0.0, 0.8}, step),
      "Boosts are not Killing in Minkowski");
    assert(GeneralRelativity.Symmetry.conformalKillingResidual(GeneralRelativity.Metrics.minkowski, GeneralRelativity.Symmetry.cartesianBoost, x, {0.0}, {0.6, 0.0, 0.8}, step) < 1.0e-10,
      "A Killing vector failed the conformal Killing test");
    // The dilation is conformal Killing but not Killing: L eta = 2 eta, div = 4.
    (L, scale) := GeneralRelativity.Symmetry.lieDerivativeMetric(GeneralRelativity.Metrics.minkowski, GeneralRelativity.Symmetry.dilation, x, {0.0}, {0.0}, step);
    assert(Tests.Assertions.maxAbsMatrix(L - 2.0 * eta) < 1.0e-10, "Lie derivative of eta along the dilation is not 2 eta");
    assert(abs(GeneralRelativity.Symmetry.divergence(GeneralRelativity.Metrics.minkowski, GeneralRelativity.Symmetry.dilation, x, {0.0}, {0.0}, step) - 4.0) < 1.0e-10,
      "Divergence of the dilation is not 4");
    assert(GeneralRelativity.Symmetry.killingResidual(GeneralRelativity.Metrics.minkowski, GeneralRelativity.Symmetry.dilation, x, {0.0}, {0.0}, step) > 0.5,
      "The dilation passed as a Killing vector");
    assert(GeneralRelativity.Symmetry.conformalKillingResidual(GeneralRelativity.Metrics.minkowski, GeneralRelativity.Symmetry.dilation, x, {0.0}, {0.0}, step) < 1.0e-10,
      "The dilation is not conformal Killing in Minkowski");
    passed := true;
  end runSymmetry;

  // ---------------------------------------------------------------------------
  function runWedge "Wedge products and interior products"
    output Boolean passed;
  protected
    constant Real rs = 2.0;
    constant Real step[4] = {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3};
    constant Real eta[4, 4] = diagonal({-1.0, 1.0, 1.0, 1.0});
    constant Real a[4] = {0.3, -1.2, 0.7, 2.0};
    constant Real b[4] = {1.1, 0.4, -0.5, 0.9};
    constant Real c[4] = {-0.6, 0.8, 1.5, 0.2};
    constant Real d[4] = {0.5, -0.3, 0.1, 1.4};
    constant Real X[4] = {0.9, -0.2, 0.4, 1.3};
    Real e0[4] = {1, 0, 0, 0};
    Real e1[4] = {0, 1, 0, 0};
    Real e2[4] = {0, 0, 1, 0};
    Real e3[4] = {0, 0, 0, 1};
    Real metrics[2, 4, 4];
    Real g[4, 4];
    Real gInv[4, 4];
    Boolean accepted;
    Real x[4];
    Real F[4, 4];
    Real F2[4, 4];
    Real H[4, 4, 4];
    Real H2[4, 4, 4];
    Real V[4, 4, 4, 4];
    Real V2[4, 4, 4, 4];
    Real eps[4, 4, 4, 4];
    Real w[4];
    Real det;
    Real err;
    Real deltaW;
    Real deltaF[4];
    Real rhs[4];
    Real T[4, 4];
    GeneralRelativity.Tensors.Curvature k;
    Real Fsq;
    Real FF[4, 4];
    Real F2Up[4, 4];
  algorithm
    // Permutation symbol and wedge algebra.
    assert(GeneralRelativity.ExteriorCalculus.leviCivitaSymbol(1, 2, 3, 4) == 1
      and GeneralRelativity.ExteriorCalculus.leviCivitaSymbol(2, 1, 3, 4) == -1
      and GeneralRelativity.ExteriorCalculus.leviCivitaSymbol(4, 1, 2, 3) == -1
      and GeneralRelativity.ExteriorCalculus.leviCivitaSymbol(3, 4, 1, 2) == 1
      and GeneralRelativity.ExteriorCalculus.leviCivitaSymbol(1, 3, 3, 4) == 0, "Permutation symbol is wrong");
    assert(Tests.Assertions.maxAbsMatrix(GeneralRelativity.ExteriorCalculus.wedge11(a, b) + GeneralRelativity.ExteriorCalculus.wedge11(b, a)) < 1.0e-15
      and Tests.Assertions.maxAbsMatrix(GeneralRelativity.ExteriorCalculus.wedge11(a, a)) < 1.0e-15, "1-form wedge is not antisymmetric");
    H := GeneralRelativity.ExteriorCalculus.wedge12(a, GeneralRelativity.ExteriorCalculus.wedge11(b, c));
    err := 0.0;
    for i in 1:4 loop
      for j in 1:4 loop
        for l in 1:4 loop
          err := max(err, max(abs(H[i, j, l] + H[j, i, l]), abs(H[i, j, l] + H[i, l, j])));
        end for;
      end for;
    end for;
    assert(err < 1.0e-14, "a ^ b ^ c is not totally antisymmetric");
    assert(maxAbs3(H - GeneralRelativity.ExteriorCalculus.wedge12(c, GeneralRelativity.ExteriorCalculus.wedge11(a, b))) < 1.0e-14
      and maxAbs3(GeneralRelativity.ExteriorCalculus.wedge12(a, GeneralRelativity.ExteriorCalculus.wedge11(a, c))) < 1.0e-14,
      "a ^ b ^ c is not cyclic or a ^ a ^ c is not zero");
    V := GeneralRelativity.ExteriorCalculus.wedge13(a, GeneralRelativity.ExteriorCalculus.wedge12(b, GeneralRelativity.ExteriorCalculus.wedge11(c, d)));
    V2 := GeneralRelativity.ExteriorCalculus.wedge22(GeneralRelativity.ExteriorCalculus.wedge11(a, b), GeneralRelativity.ExteriorCalculus.wedge11(c, d));
    assert(maxAbs4(V - V2) < 1.0e-13, "Wedge product is not associative: a ^ (b ^ c ^ d) differs from (a ^ b) ^ (c ^ d)");
    det := 0.0;
    for i in 1:4 loop
      for j in 1:4 loop
        for l in 1:4 loop
          for m in 1:4 loop
            det := det + GeneralRelativity.ExteriorCalculus.leviCivitaSymbol(i, j, l, m) * a[i] * b[j] * c[l] * d[m];
          end for;
        end for;
      end for;
    end for;
    assert(abs(det) > 0.1 and maxAbs4(V - det * GeneralRelativity.ExteriorCalculus.leviCivitaTensor(eta)) < 1.0e-13,
      "a ^ b ^ c ^ d is not det[a; b; c; d] times the permutation symbol");
    F := GeneralRelativity.ExteriorCalculus.wedge11(a, b) + GeneralRelativity.ExteriorCalculus.wedge11(c, d);
    F2 := GeneralRelativity.ExteriorCalculus.wedge11(b, d) - 0.5 * GeneralRelativity.ExteriorCalculus.wedge11(a, c);
    assert(maxAbs4(GeneralRelativity.ExteriorCalculus.wedge22(F, F2) - GeneralRelativity.ExteriorCalculus.wedge22(F2, F)) < 1.0e-13,
      "2-form wedge is not symmetric");

    // Interior product is an antiderivation and nilpotent.
    assert(Tests.Assertions.maxAbsVector(GeneralRelativity.ExteriorCalculus.interior2(X, GeneralRelativity.ExteriorCalculus.wedge11(a, b))
      - (GeneralRelativity.ExteriorCalculus.interior1(X, a) * b - GeneralRelativity.ExteriorCalculus.interior1(X, b) * a)) < 1.0e-14,
      "i_X (a ^ b) is not (i_X a) b - (i_X b) a");
    assert(Tests.Assertions.maxAbsMatrix(GeneralRelativity.ExteriorCalculus.interior3(X, GeneralRelativity.ExteriorCalculus.wedge12(a, F))
      - (GeneralRelativity.ExteriorCalculus.interior1(X, a) * F - GeneralRelativity.ExteriorCalculus.wedge11(a, GeneralRelativity.ExteriorCalculus.interior2(X, F)))) < 1.0e-13,
      "i_X (a ^ F) is not (i_X a) F - a ^ i_X F");
    assert(Tests.Assertions.maxAbsVector(GeneralRelativity.ExteriorCalculus.interior2(X, GeneralRelativity.ExteriorCalculus.interior3(X, H))) < 1.0e-13
      and Tests.Assertions.maxAbsMatrix(GeneralRelativity.ExteriorCalculus.interior3(X, GeneralRelativity.ExteriorCalculus.interior4(X, V))) < 1.0e-13,
      "i_X i_X is not zero");

    passed := true;
  end runWedge;

  function runHodge "Levi-Civita tensor and Hodge duals in Minkowski and Schwarzschild"
    output Boolean passed;
  protected
    constant Real rs = 2.0;
    constant Real step[4] = {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3};
    constant Real eta[4, 4] = diagonal({-1.0, 1.0, 1.0, 1.0});
    constant Real a[4] = {0.3, -1.2, 0.7, 2.0};
    constant Real b[4] = {1.1, 0.4, -0.5, 0.9};
    constant Real c[4] = {-0.6, 0.8, 1.5, 0.2};
    constant Real d[4] = {0.5, -0.3, 0.1, 1.4};
    constant Real X[4] = {0.9, -0.2, 0.4, 1.3};
    Real e0[4] = {1, 0, 0, 0};
    Real e1[4] = {0, 1, 0, 0};
    Real e2[4] = {0, 0, 1, 0};
    Real e3[4] = {0, 0, 0, 1};
    Real metrics[2, 4, 4];
    Real g[4, 4];
    Real gInv[4, 4];
    Boolean accepted;
    Real x[4];
    Real F[4, 4];
    Real F2[4, 4];
    Real H[4, 4, 4];
    Real H2[4, 4, 4];
    Real V[4, 4, 4, 4];
    Real V2[4, 4, 4, 4];
    Real eps[4, 4, 4, 4];
    Real w[4];
    Real det;
    Real err;
    Real deltaW;
    Real deltaF[4];
    Real rhs[4];
    Real T[4, 4];
    GeneralRelativity.Tensors.Curvature k;
    Real Fsq;
    Real FF[4, 4];
    Real F2Up[4, 4];
  algorithm
    F := GeneralRelativity.ExteriorCalculus.wedge11(a, b) + GeneralRelativity.ExteriorCalculus.wedge11(c, d);
    F2 := GeneralRelativity.ExteriorCalculus.wedge11(b, d) - 0.5 * GeneralRelativity.ExteriorCalculus.wedge11(a, c);
    H := GeneralRelativity.ExteriorCalculus.wedge12(a, GeneralRelativity.ExteriorCalculus.wedge11(b, c));
    V := GeneralRelativity.ExteriorCalculus.wedge13(a, GeneralRelativity.ExteriorCalculus.wedge12(b, GeneralRelativity.ExteriorCalculus.wedge11(c, d)));
    // Hodge duals in Minkowski: *(dt ^ dx) = -dy ^ dz, *dt = -dx ^ dy ^ dz.
    assert(Tests.Assertions.maxAbsMatrix(GeneralRelativity.ExteriorCalculus.hodge2(GeneralRelativity.ExteriorCalculus.wedge11(e0, e1), eta)
      + GeneralRelativity.ExteriorCalculus.wedge11(e2, e3)) < 1.0e-15, "*(dt ^ dx) is not -dy ^ dz");
    assert(maxAbs3(GeneralRelativity.ExteriorCalculus.hodge1(e0, eta)
      + GeneralRelativity.ExteriorCalculus.wedge12(e1, GeneralRelativity.ExteriorCalculus.wedge11(e2, e3))) < 1.0e-15, "*dt is not -dx ^ dy ^ dz");

    // ** = (-1)^(p(4-p)+1), the Hodge of the volume form and alpha ^ *beta = <alpha, beta> eps,
    // in Minkowski and in Schwarzschild spherical coordinates.
    x := {0.3, 5.0 * rs, 1.1, 0.4};
    metrics[1, :, :] := eta;
    metrics[2, :, :] := GeneralRelativity.Metrics.schwarzschild(x, {rs});
    for m in 1:2 loop
      g := metrics[m, :, :];
      (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(g);
      eps := GeneralRelativity.ExteriorCalculus.leviCivitaTensor(g);
      F2Up := gInv * F2 * gInv;
      assert(abs(GeneralRelativity.ExteriorCalculus.hodge4(eps, g) + 1.0) < 1.0e-12, "*eps is not -1");
      assert(maxAbs4(GeneralRelativity.ExteriorCalculus.hodge0(1.0, g) - eps) < 1.0e-15, "*1 is not the volume form");
      assert(abs(GeneralRelativity.ExteriorCalculus.hodge4(GeneralRelativity.ExteriorCalculus.hodge0(2.5, g), g) + 2.5) < 1.0e-12, "** is not -1 on 0-forms");
      assert(Tests.Assertions.maxAbsVector(GeneralRelativity.ExteriorCalculus.hodge3(GeneralRelativity.ExteriorCalculus.hodge1(a, g), g) - a) < 1.0e-12,
        "** is not +1 on 1-forms");
      assert(Tests.Assertions.maxAbsMatrix(GeneralRelativity.ExteriorCalculus.hodge2(GeneralRelativity.ExteriorCalculus.hodge2(F, g), g) + F) < 1.0e-12,
        "** is not -1 on 2-forms");
      assert(maxAbs3(GeneralRelativity.ExteriorCalculus.hodge1(GeneralRelativity.ExteriorCalculus.hodge3(H, g), g) - H) < 1.0e-12,
        "** is not +1 on 3-forms");
      assert(maxAbs4(GeneralRelativity.ExteriorCalculus.hodge0(GeneralRelativity.ExteriorCalculus.hodge4(V, g), g) + V) < 1.0e-11,
        "** is not -1 on 4-forms");
      assert(maxAbs4(GeneralRelativity.ExteriorCalculus.wedge13(a, GeneralRelativity.ExteriorCalculus.hodge1(b, g)) - (a * gInv * b) * eps) < 1.0e-11,
        "a ^ *b is not <a, b> eps");
      assert(maxAbs4(GeneralRelativity.ExteriorCalculus.wedge22(F, GeneralRelativity.ExteriorCalculus.hodge2(F2, g))
        - 0.5 * sum(F[i, j] * F2Up[i, j] for i in 1:4, j in 1:4) * eps) < 1.0e-11,
        "F ^ *G is not (1/2) F_ab G^ab eps");
      assert(maxAbs3(GeneralRelativity.ExteriorCalculus.interior4(X, eps) - GeneralRelativity.ExteriorCalculus.hodge1(g * X, g)) < 1.0e-11,
        "i_X eps is not the dual of the lowered vector");
    end for;
    assert(abs(GeneralRelativity.ExteriorCalculus.volumeElement(metrics[2, :, :]) - x[2]^2 * sin(x[3])) < 1.0e-12,
      "Schwarzschild volume element is not r^2 sin(theta)");

    passed := true;
  end runHodge;

  function runDerivatives "Exterior derivative and codifferential"
    output Boolean passed;
  protected
    constant Real rs = 2.0;
    constant Real step[4] = {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3};
    constant Real eta[4, 4] = diagonal({-1.0, 1.0, 1.0, 1.0});
    constant Real a[4] = {0.3, -1.2, 0.7, 2.0};
    constant Real b[4] = {1.1, 0.4, -0.5, 0.9};
    constant Real c[4] = {-0.6, 0.8, 1.5, 0.2};
    constant Real d[4] = {0.5, -0.3, 0.1, 1.4};
    constant Real X[4] = {0.9, -0.2, 0.4, 1.3};
    Real e0[4] = {1, 0, 0, 0};
    Real e1[4] = {0, 1, 0, 0};
    Real e2[4] = {0, 0, 1, 0};
    Real e3[4] = {0, 0, 0, 1};
    Real metrics[2, 4, 4];
    Real g[4, 4];
    Real gInv[4, 4];
    Boolean accepted;
    Real x[4];
    Real F[4, 4];
    Real F2[4, 4];
    Real H[4, 4, 4];
    Real H2[4, 4, 4];
    Real V[4, 4, 4, 4];
    Real V2[4, 4, 4, 4];
    Real eps[4, 4, 4, 4];
    Real w[4];
    Real det;
    Real err;
    Real deltaW;
    Real deltaF[4];
    Real rhs[4];
    Real T[4, 4];
    GeneralRelativity.Tensors.Curvature k;
    Real Fsq;
    Real FF[4, 4];
    Real F2Up[4, 4];
  algorithm
    // Exterior derivative: closed form for testOneForm, and d d = 0 on 0-, 1- and 2-forms.
    x := {0.4, 1.2, -0.7, 0.9};
    F := GeneralRelativity.ExteriorCalculus.exteriorDerivative1(testOneForm, x, {0.0}, step);
    assert(abs(F[1, 2] - (2.0 * x[1] - x[3] * cos(x[2]))) < 1.0e-6 and abs(F[2, 4] - sin(x[4])) < 1.0e-6
      and abs(F[3, 4] - (exp(0.1 * x[1]) - x[2])) < 1.0e-6 and abs(F[2, 1] + F[1, 2]) < 1.0e-15 and abs(F[1, 4] - 0.1 * exp(0.1 * x[1]) * x[3]) < 1.0e-6,
      "Exterior derivative of the test 1-form is wrong");
    w := GeneralRelativity.ExteriorCalculus.exteriorDerivative0(testScalar, x, {0.0}, step);
    assert(abs(w[1] - x[2] * cos(x[1] * x[2])) < 1.0e-6, "Gradient of the test scalar is wrong");
    assert(Tests.Assertions.maxAbsMatrix(GeneralRelativity.ExteriorCalculus.exteriorDerivative1(gradientOfTestScalar, x, {0.0}, step)) < 1.0e-8,
      "d d f is not zero");
    assert(maxAbs3(GeneralRelativity.ExteriorCalculus.exteriorDerivative2(dOfTestOneForm, x, {0.0}, step)) < 1.0e-8,
      "d d A is not zero");
    assert(maxAbs4(GeneralRelativity.ExteriorCalculus.exteriorDerivative3(dOfTestTwoForm, x, {0.0}, step)) < 1.0e-8,
      "d d F is not zero");

    // Codifferential: delta = -div; delta(x_a dx^a) = -4 in Minkowski; delta = *d* on Schwarzschild.
    assert(abs(GeneralRelativity.ExteriorCalculus.codifferential1(GeneralRelativity.Metrics.minkowski, positionOneForm, x, {0.0}, {0.0}, step) + 4.0) < 1.0e-10,
      "Codifferential of the position 1-form is not -4");
    x := {0.3, 5.0 * rs, 1.1, 0.4};
    g := GeneralRelativity.Metrics.schwarzschild(x, {rs});
    deltaW := GeneralRelativity.ExteriorCalculus.codifferential1(GeneralRelativity.Metrics.schwarzschild, testOneForm, x, {rs}, {rs}, step);
    assert(abs(deltaW) > 0.01 and abs(deltaW - GeneralRelativity.ExteriorCalculus.hodge4(GeneralRelativity.ExteriorCalculus.exteriorDerivative3(dualOfTestOneForm, x, {rs}, step), g)) < 1.0e-6,
      "delta w differs from *d*w");
    deltaF := GeneralRelativity.ExteriorCalculus.codifferential2(GeneralRelativity.Metrics.schwarzschild, testTwoForm, x, {rs}, {rs}, step);
    rhs := GeneralRelativity.ExteriorCalculus.hodge3(GeneralRelativity.ExteriorCalculus.exteriorDerivative2(dualOfTestTwoForm, x, {rs}, step), g);
    assert(Tests.Assertions.maxAbsVector(deltaF) > 0.01 and Tests.Assertions.maxAbsVector(deltaF - rhs) < 1.0e-6,
      "delta F differs from *d*F");

    passed := true;
  end runDerivatives;

  function runMaxwell "Coulomb field on Reissner-Nordstrom"
    output Boolean passed;
  protected
    constant Real rs = 2.0;
    constant Real step[4] = {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3};
    constant Real eta[4, 4] = diagonal({-1.0, 1.0, 1.0, 1.0});
    constant Real a[4] = {0.3, -1.2, 0.7, 2.0};
    constant Real b[4] = {1.1, 0.4, -0.5, 0.9};
    constant Real c[4] = {-0.6, 0.8, 1.5, 0.2};
    constant Real d[4] = {0.5, -0.3, 0.1, 1.4};
    constant Real X[4] = {0.9, -0.2, 0.4, 1.3};
    Real e0[4] = {1, 0, 0, 0};
    Real e1[4] = {0, 1, 0, 0};
    Real e2[4] = {0, 0, 1, 0};
    Real e3[4] = {0, 0, 0, 1};
    Real metrics[2, 4, 4];
    Real g[4, 4];
    Real gInv[4, 4];
    Boolean accepted;
    Real x[4];
    Real F[4, 4];
    Real F2[4, 4];
    Real H[4, 4, 4];
    Real H2[4, 4, 4];
    Real V[4, 4, 4, 4];
    Real V2[4, 4, 4, 4];
    Real eps[4, 4, 4, 4];
    Real w[4];
    Real det;
    Real err;
    Real deltaW;
    Real deltaF[4];
    Real rhs[4];
    Real T[4, 4];
    GeneralRelativity.Tensors.Curvature k;
    Real Fsq;
    Real FF[4, 4];
    Real F2Up[4, 4];
  algorithm
    // Maxwell on Reissner-Nordstrom: the Coulomb field (Q/r^2) dt ^ dr is closed and co-closed,
    // and with rQ = Q (geometrized Gaussian units) it sources the metric: G_ab = 2 (F_ac F_b^c - g_ab F^2 / 4).
    x := {0.1, 6.0, 1.0, 0.3};
    assert(maxAbs3(GeneralRelativity.ExteriorCalculus.exteriorDerivative2(radialFieldStrength, x, {0.8, 2.0}, step)) < 1.0e-12,
      "Coulomb field is not closed");
    deltaF := GeneralRelativity.ExteriorCalculus.codifferential2(reissnerNordstrom, radialFieldStrength, x, {rs, 0.8}, {0.8, 2.0}, step);
    assert(Tests.Assertions.maxAbsVector(deltaF) < 1.0e-8, "Coulomb field is not co-closed on Reissner-Nordstrom");
    deltaF := GeneralRelativity.ExteriorCalculus.codifferential2(reissnerNordstrom, radialFieldStrength, x, {rs, 0.8}, {0.8, 3.0}, step);
    // (Q/r^3) dx0 ^ dr is not co-closed: (delta F)_0 = -f Q / r^4 with f = 1 - rs/r + rQ^2/r^2.
    assert(abs(deltaF[1] + (1.0 - rs / x[2] + 0.8^2 / x[2]^2) * 0.8 / x[2]^4) < 1.0e-9 and abs(deltaF[1]) > 1.0e-4,
      "Codifferential of a 1/r^3 radial field is not -f Q / r^4: " + String(deltaF[1]));
    k := GeneralRelativity.Metrics.curvatureAt(reissnerNordstrom, x, {rs, 0.8}, step);
    F := radialFieldStrength(x, {0.8, 2.0});
    FF := F * k.gInv * transpose(F);
    F2Up := k.gInv * F * k.gInv;
    Fsq := sum(F[i, j] * F2Up[i, j] for i in 1:4, j in 1:4);
    T := 2.0 * (FF - 0.25 * Fsq * k.g);
    assert(Tests.Assertions.maxAbsMatrix(k.einstein - T) < 1.0e-6 and Tests.Assertions.maxAbsMatrix(T) > 1.0e-3,
      "Reissner-Nordstrom Einstein tensor is not the Maxwell stress-energy of the Coulomb field");
    passed := true;
  end runMaxwell;

  // ---------------------------------------------------------------------------
  function runCartan "Spin connection and Cartan's structure equations on Schwarzschild"
    output Boolean passed;
  protected
    constant Real rs = 2.0;
    constant Real step[4] = {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3};
    constant Real eta[4] = {-1.0, 1.0, 1.0, 1.0};
    Real x[4] = {0.3, 5.0 * rs, 1.1, 0.4};
    Real g[4, 4];
    Real e[4, 4];
    Real theta[4, 4];
    Real omega[4, 4, 4];
    Real T[4, 4, 4];
    Real residual[4, 4, 4, 4];
    Real relative;
    Real err;
    Real f;
    Real r;
  algorithm
    r := x[2];
    f := 1.0 - rs / r;
    g := GeneralRelativity.Metrics.schwarzschild(x, {rs});
    e := GeneralRelativity.ExteriorCalculus.schwarzschildTetrad(x, {rs});
    assert(GeneralRelativity.ExteriorCalculus.orthonormalityResidual(g, e) < 1.0e-14, "Schwarzschild tetrad is not orthonormal");
    theta := GeneralRelativity.ExteriorCalculus.coframe(g, e);
    assert(Tests.Assertions.maxAbsMatrix(theta * transpose(e) - identity(4)) < 1.0e-14, "Coframe is not dual to the tetrad");
    assert(Tests.Assertions.maxAbsMatrix(transpose(theta) * diagonal(eta) * theta - g) < 1.0e-12, "Coframe does not rebuild the metric");

    // Closed-form connection forms of the static frame.
    omega := GeneralRelativity.ExteriorCalculus.spinConnection(GeneralRelativity.Metrics.schwarzschild, GeneralRelativity.ExteriorCalculus.schwarzschildTetrad, x, {rs}, {rs}, step);
    assert(abs(omega[1, 2, 1] - rs / (2.0 * r^2)) < 1.0e-8 and abs(omega[2, 1, 1] - rs / (2.0 * r^2)) < 1.0e-8,
      "omega^0_1 is not (f'/2) dx0");
    assert(abs(omega[3, 2, 3] - sqrt(f)) < 1.0e-6 and abs(omega[4, 2, 4] - sqrt(f) * sin(x[3])) < 1.0e-6
      and abs(omega[4, 3, 4] - cos(x[3])) < 1.0e-6, "Angular connection forms of the static frame are wrong: "
      + String(omega[3, 2, 3] - sqrt(f)) + " " + String(omega[4, 2, 4] - sqrt(f) * sin(x[3])) + " " + String(omega[4, 3, 4] - cos(x[3])));
    err := 0.0;
    for a in 1:4 loop
      for b in 1:4 loop
        for mu in 1:4 loop
          err := max(err, abs(eta[a] * omega[a, b, mu] + eta[b] * omega[b, a, mu]));
        end for;
      end for;
    end for;
    // Truncation error of d_theta e is about 1e-7 relative; the coframe theta^3 = r sin(theta) dphi multiplies it by r.
    assert(err < 1.0e-6 * r, "omega_ab is not antisymmetric: " + String(err));

    // First and second structure equations for the static frame and a locally boosted and rotated frame.
    // Both residuals are O(step^2) truncation error (halving the step divides them by 4); the torsion
    // components carry the size r of the coframe, the curvature residual is relative.
    T := GeneralRelativity.ExteriorCalculus.cartanTorsion(GeneralRelativity.Metrics.schwarzschild, GeneralRelativity.ExteriorCalculus.schwarzschildTetrad, x, {rs}, {rs}, step);
    assert(maxAbs3(T) < 2.0e-6 * r, "d theta + omega ^ theta is not zero for the static frame: " + String(maxAbs3(T)));
    (residual, relative) := GeneralRelativity.ExteriorCalculus.cartanCurvatureResidual(GeneralRelativity.Metrics.schwarzschild, GeneralRelativity.ExteriorCalculus.schwarzschildTetrad, x, {rs}, {rs}, step);
    assert(relative < 1.0e-5, "Cartan curvature 2-forms differ from the Riemann tensor for the static frame: " + String(relative));

    e := rotatedTetrad(x, {rs});
    assert(GeneralRelativity.ExteriorCalculus.orthonormalityResidual(g, e) < 1.0e-12, "Rotated tetrad is not orthonormal");
    omega := GeneralRelativity.ExteriorCalculus.spinConnection(GeneralRelativity.Metrics.schwarzschild, rotatedTetrad, x, {rs}, {rs}, step);
    err := 0.0;
    for a in 1:4 loop
      for b in 1:4 loop
        for mu in 1:4 loop
          err := max(err, abs(eta[a] * omega[a, b, mu] + eta[b] * omega[b, a, mu]));
        end for;
      end for;
    end for;
    assert(err < 1.0e-6 * r, "omega_ab of the rotated frame is not antisymmetric: " + String(err));
    T := GeneralRelativity.ExteriorCalculus.cartanTorsion(GeneralRelativity.Metrics.schwarzschild, rotatedTetrad, x, {rs}, {rs}, step);
    assert(maxAbs3(T) < 2.0e-6 * r, "d theta + omega ^ theta is not zero for the rotated frame: " + String(maxAbs3(T)));
    (residual, relative) := GeneralRelativity.ExteriorCalculus.cartanCurvatureResidual(GeneralRelativity.Metrics.schwarzschild, rotatedTetrad, x, {rs}, {rs}, step);
    assert(relative < 1.0e-5, "Cartan curvature 2-forms differ from the Riemann tensor for the rotated frame: " + String(relative));
    passed := true;
  end runCartan;

  parameter Boolean foliationPassed = runFoliation();
  parameter Boolean symmetryPassed = runSymmetry();
  parameter Boolean wedgePassed = runWedge();
  parameter Boolean hodgePassed = runHodge();
  parameter Boolean derivativesPassed = runDerivatives();
  parameter Boolean maxwellPassed = runMaxwell();
  parameter Boolean cartanPassed = runCartan();
equation
  assert(foliationPassed and symmetryPassed and wedgePassed and hodgePassed
    and derivativesPassed and maxwellPassed and cartanPassed, "GeneralRelativity geometry assertions did not complete");
  annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 1.0e-8, Interval = 0.1));
end GeneralRelativityGeometryTests;
