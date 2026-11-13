within Tests;
model GeneralRelativityCurvatureTests "Weyl decomposition, curvature invariants, Cotton tensor, tetrads, Newman-Penrose scalars, Petrov types and energy conditions"
  function staticSphere "Non-vacuum static spherically symmetric test metric in {ct, r, theta, phi}"
    extends GeneralRelativity.Interfaces.MetricFunction;
  algorithm
    g := zeros(4, 4);
    g[1, 1] := -(1.0 + x[2]^2);
    g[2, 2] := 1.0 / (1.0 - x[2]^2 / 4.0);
    g[3, 3] := x[2]^2;
    g[4, 4] := x[2]^2 * sin(x[3])^2;
  end staticSphere;

  function conformallyFlat "exp(2 sigma(x)) times the Minkowski metric in Cartesian coordinates"
    extends GeneralRelativity.Interfaces.MetricFunction;
  protected
    Real sigma = 0.1 * x[1] + 0.05 * x[2]^2 - 0.08 * x[3] * x[4] + 0.03 * x[4];
  algorithm
    g := exp(2.0 * sigma) * diagonal({-1.0, 1.0, 1.0, 1.0});
  end conformallyFlat;

  function schwarzschildDeSitter "Schwarzschild-de Sitter metric in {ct, r, theta, phi}; p = {rs, Lambda}"
    extends GeneralRelativity.Interfaces.MetricFunction;
  protected
    Real f = 1.0 - p[1] / x[2] - p[2] * x[2]^2 / 3.0;
  algorithm
    g := zeros(4, 4);
    g[1, 1] := -f;
    g[2, 2] := 1.0 / f;
    g[3, 3] := x[2]^2;
    g[4, 4] := x[2]^2 * sin(x[3])^2;
  end schwarzschildDeSitter;

  function ppWave "Vacuum plane wave -du dv + H du^2 + dx^2 + dy^2, u = t - z, v = t + z, H = (x^2 - y^2)(1 + sin(u)/2), in {ct, x, y, z}"
    extends GeneralRelativity.Interfaces.MetricFunction;
  protected
    Real H = (x[2]^2 - x[3]^2) * (1.0 + 0.5 * sin(x[1] - x[4]));
  algorithm
    g := diagonal({-1.0, 1.0, 1.0, 1.0});
    g[1, 1] := g[1, 1] + H;
    g[1, 4] := -H;
    g[4, 1] := -H;
    g[4, 4] := g[4, 4] + H;
  end ppWave;

  function genericMetric "Minkowski plus a quadratic perturbation without symmetries, for an algebraically general Weyl tensor"
    extends GeneralRelativity.Interfaces.MetricFunction;
  algorithm
    g := diagonal({-1.0, 1.0, 1.0, 1.0});
    for a in 1:4 loop
      for b in 1:4 loop
        g[a, b] := g[a, b] + 0.1 * sum(sin(1.3 * (a + b) + 0.7 * a * b + 1.1 * c + 0.5 * c * d + 0.9 * d) * x[c] * x[d] for c in 1:4, d in 1:4);
      end for;
    end for;
  end genericMetric;

  function kulkarniNomizu "Algebraic curvature tensor (g wedge S)_abcd = g_ac S_bd - g_ad S_bc - g_bc S_ad + g_bd S_ac"
    input Real g[:, size(g, 1)];
    input Real S[size(g, 1), size(g, 1)];
    output Real R[size(g, 1), size(g, 1), size(g, 1), size(g, 1)];
  algorithm
    for a in 1:size(g, 1) loop
      for b in 1:size(g, 1) loop
        for c in 1:size(g, 1) loop
          for d in 1:size(g, 1) loop
            R[a, b, c, d] := g[a, c] * S[b, d] - g[a, d] * S[b, c] - g[b, c] * S[a, d] + g[b, d] * S[a, c];
          end for;
        end for;
      end for;
    end for;
  end kulkarniNomizu;

  function maxAbs4 "Maximum absolute element of a rank-4 array"
    input Real A[:, :, :, :];
    output Real m;
  algorithm
    m := max(abs(A[a, b, c, d]) for a in 1:size(A, 1), b in 1:size(A, 2), c in 1:size(A, 3), d in 1:size(A, 4));
  end maxAbs4;

  function maxAbs3 "Maximum absolute element of a rank-3 array"
    input Real A[:, :, :];
    output Real m;
  algorithm
    m := max(abs(A[a, b, c]) for a in 1:size(A, 1), b in 1:size(A, 2), c in 1:size(A, 3));
  end maxAbs3;

  function weylTraceError "Largest trace g^ac C_abcd of a rank-4 tensor"
    input Real gInv[4, 4];
    input Real C[4, 4, 4, 4];
    output Real m;
  algorithm
    m := max(abs(sum(gInv[a, c] * C[a, b, c, d] for a in 1:4, c in 1:4)) for b in 1:4, d in 1:4);
  end weylTraceError;

  function tetradError "Largest deviation of e g e^T from diag(-1, 1, 1, 1)"
    input Real g[4, 4];
    input Real e[4, 4];
    output Real m;
  algorithm
    m := Tests.Assertions.maxAbsMatrix(e * g * transpose(e) - diagonal({-1.0, 1.0, 1.0, 1.0}));
  end tetradError;

  function run
    output Boolean passed;
  protected
    constant Real pi = 3.141592653589793;
    constant Real rs = 2.0;
    constant Real Lambda = 0.01;
    constant Real step[4] = {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3};
    constant Real eta[4, 4] = [-1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
    GeneralRelativity.Tensors.Curvature bh;
    GeneralRelativity.Tensors.Curvature bh2;
    GeneralRelativity.Tensors.Curvature kerr;
    GeneralRelativity.Tensors.Curvature star;
    GeneralRelativity.Tensors.Curvature sds;
    GeneralRelativity.Tensors.Curvature conf;
    GeneralRelativity.Tensors.Curvature wave;
    GeneralRelativity.Tensors.Curvature generic;
    GeneralRelativity.Tensors.Invariants inv;
    GeneralRelativity.Tensors.Invariants inv2;
    Real x[4];
    Real r;
    Real m;
    Real C[4, 4, 4, 4];
    Real F[4, 4, 4, 4];
    Real P[4, 4];
    Real cotton[4, 4, 4];
    Real cottonTrace;
    Real e[4, 4];
    Real g[4, 4];
    Real gInv[4, 4];
    Boolean accepted;
    Real l[4];
    Real n[4];
    Real mRe[4];
    Real mIm[4];
    Real psiRe[5];
    Real psiIm[5];
    Real IRe;
    Real IIm;
    Real JRe;
    Real JIm;
    Real I3Re;
    Real I3Im;
    Real J2Re;
    Real J2Im;
    Real E3[3, 3];
    Real B3[3, 3];
    Real lambda;
    Real cubic[6];
    Real expected;
    Real u[4];
    Real T[4, 4];
    Real k[4];
    Boolean nec;
    Boolean wec;
    Boolean sec;
    Boolean dec;
    Real rho;
    Real pr[3];
    Boolean typeI;
    Real g3[3, 3];
    Real S3[3, 3];
    Real R3[3, 3, 3, 3];
    Real g3Inv[3, 3];
    Real Ric3[3, 3];
    Real S4[4, 4];
    Real R4[4, 4, 4, 4];
    Real Ric4[4, 4];
  algorithm
    // --- Algebraic decomposition: a pure g wedge S curvature has Schouten tensor S and no Weyl part. ---
    g3 := [2.0, 0.3, -0.1; 0.3, 1.5, 0.2; -0.1, 0.2, 1.0];
    S3 := [0.4, -0.2, 0.7; -0.2, 1.1, 0.05; 0.7, 0.05, -0.6];
    R3 := kulkarniNomizu(g3, S3);
    (g3Inv, accepted) := GeneralRelativity.Tensors.inverseMetric(g3);
    Ric3 := GeneralRelativity.Tensors.ricci(g3Inv, R3);
    assert(Tests.Assertions.maxAbsMatrix(GeneralRelativity.Tensors.schouten(g3, Ric3, GeneralRelativity.Tensors.ricciScalar(g3Inv, Ric3)) - S3) < 1.0e-12,
      "3D Schouten tensor of g wedge S is not S");
    assert(maxAbs4(GeneralRelativity.Tensors.weyl(g3, R3, Ric3, GeneralRelativity.Tensors.ricciScalar(g3Inv, Ric3))) < 1.0e-12,
      "Weyl tensor is not identically zero in three dimensions");
    g := [-1.2, 0.1, 0.0, 0.2; 0.1, 1.3, 0.1, 0.0; 0.0, 0.1, 0.9, -0.1; 0.2, 0.0, -0.1, 1.1];
    S4 := [0.3, 0.1, -0.2, 0.0; 0.1, -0.5, 0.4, 0.3; -0.2, 0.4, 0.8, -0.1; 0.0, 0.3, -0.1, 0.2];
    R4 := kulkarniNomizu(g, S4);
    (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(g);
    Ric4 := GeneralRelativity.Tensors.ricci(gInv, R4);
    assert(Tests.Assertions.maxAbsMatrix(GeneralRelativity.Tensors.schouten(g, Ric4, GeneralRelativity.Tensors.ricciScalar(gInv, Ric4)) - S4) < 1.0e-12,
      "4D Schouten tensor of g wedge S is not S");
    assert(maxAbs4(GeneralRelativity.Tensors.weyl(g, R4, Ric4, GeneralRelativity.Tensors.ricciScalar(gInv, Ric4))) < 1.0e-12,
      "Weyl tensor of a pure g wedge S curvature is not zero");
    assert(abs(GeneralRelativity.Tensors.ricciScalar(gInv, GeneralRelativity.Tensors.traceFreeRicci(g, Ric4, GeneralRelativity.Tensors.ricciScalar(gInv, Ric4)))) < 1.0e-12,
      "traceFreeRicci is not traceless");

    // --- Schwarzschild: Weyl = Riemann, invariants, cubic invariant, frame components, Psi2. ---
    r := 5.0 * rs;
    x := {0.0, r, 1.1, 0.4};
    bh := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.schwarzschild, x, {rs}, step);
    C := GeneralRelativity.Tensors.weyl(bh.g, bh.riemann, bh.ricci, bh.ricciScalar);
    // Compare in the orthonormal frame, where all components share one scale.
    e := GeneralRelativity.Frames.orthonormalTetrad(bh.g);
    assert(maxAbs4(GeneralRelativity.Frames.frameComponents4(C - bh.riemann, e)) < 1.0e-5 * maxAbs4(GeneralRelativity.Frames.frameComponents4(bh.riemann, e)),
      "Schwarzschild Weyl tensor differs from Riemann");
    assert(weylTraceError(bh.gInv, C) < 1.0e-6 * maxAbs4(C), "Schwarzschild Weyl tensor is not traceless");
    assert(abs(GeneralRelativity.Tensors.fullContraction(bh.gInv, bh.riemann, bh.riemann) / bh.kretschmann - 1.0) < 1.0e-12,
      "fullContraction(R, R) is not the Kretschmann scalar");
    inv := GeneralRelativity.Tensors.curvatureInvariants(bh);
    assert(abs(inv.weylSquared / bh.kretschmann - 1.0) < 1.0e-6, "Schwarzschild Weyl^2 is not the Kretschmann scalar");
    assert(abs(inv.gaussBonnet / bh.kretschmann - 1.0) < 1.0e-6, "Schwarzschild Gauss-Bonnet density is not the Kretschmann scalar");
    assert(abs(inv.ricciSquared) < 1.0e-6 * bh.kretschmann, "Schwarzschild R_ab R^ab is not zero");
    // Cubic invariant: Weyl bivector eigenvalues are 2M/r^3 (x2) and -M/r^3 (x4) in R_ab^ab, so tr M^3 = 96 M^3/r^9.
    expected := 12.0 * rs^3 / r^9;
    assert(abs(inv.cubic[6] / expected - 1.0) < 1.0e-5, "Schwarzschild R_ab^cd R_cd^ef R_ef^ab is not 12 rs^3 / r^9");
    assert(abs(inv.weylCubed / expected - 1.0) < 1.0e-5, "Schwarzschild Weyl cubed is not 12 rs^3 / r^9");
    // Ricci-containing invariants are finite-difference noise (relative 1e-5) times curvature squared.
    assert(Tests.Assertions.maxAbsVector(inv.cubic[1:5]) < 1.0e-5 * bh.kretschmann^1.5, "Schwarzschild Ricci cubic invariants are not zero");
    bh2 := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.schwarzschild, {0.0, 3.0 * rs, 0.7, 0.0}, {rs}, step);
    inv2 := GeneralRelativity.Tensors.curvatureInvariants(bh2);
    assert(abs(inv2.cubic[6] / inv.cubic[6] - (5.0 / 3.0)^9) < 1.0e-5 * (5.0 / 3.0)^9, "Schwarzschild cubic invariant does not scale as r^-9");
    // Static tetrad: radial tidal component R_0101 = -rs / r^3.
    e := GeneralRelativity.Frames.orthonormalTetrad(bh.g);
    assert(tetradError(bh.g, e) < 1.0e-12, "Schwarzschild tetrad is not orthonormal");
    F := GeneralRelativity.Frames.frameComponents4(bh.riemann, e);
    assert(abs(F[1, 2, 1, 2] / (-rs / r^3) - 1.0) < 1.0e-5, "Schwarzschild frame component R_0101 is not -rs/r^3");
    assert(abs(F[1, 3, 1, 3] / (rs / (2.0 * r^3)) - 1.0) < 1.0e-5, "Schwarzschild frame component R_0202 is not rs/(2 r^3)");
    assert(Tests.Assertions.maxAbsMatrix(GeneralRelativity.Frames.frameComponents2(bh.g, e) - eta) < 1.0e-12, "Frame metric is not eta");
    (l, n, mRe, mIm) := GeneralRelativity.Frames.nullTetrad(e);
    assert(abs(l * bh.g * l) < 1.0e-12 and abs(n * bh.g * n) < 1.0e-12 and abs(l * bh.g * n + 1.0) < 1.0e-12,
      "l and n are not null with l.n = -1");
    assert(abs(mRe * bh.g * mRe - mIm * bh.g * mIm) < 1.0e-12 and abs(mRe * bh.g * mIm) < 1.0e-12
      and abs(mRe * bh.g * mRe + mIm * bh.g * mIm - 1.0) < 1.0e-12, "m is not null with m.mbar = 1");
    assert(abs(l * bh.g * mRe) < 1.0e-12 and abs(l * bh.g * mIm) < 1.0e-12 and abs(n * bh.g * mRe) < 1.0e-12 and abs(n * bh.g * mIm) < 1.0e-12,
      "m is not orthogonal to l and n");
    (psiRe, psiIm) := GeneralRelativity.Frames.weylScalars(C, l, n, mRe, mIm);
    assert(abs(psiRe[3] / (-rs / (2.0 * r^3)) - 1.0) < 1.0e-5 and abs(psiIm[3]) < 1.0e-6 * abs(psiRe[3]),
      "Schwarzschild Psi2 is not -M/r^3");
    assert(max(abs(psiRe[{1, 2, 4, 5}])) < 1.0e-5 * abs(psiRe[3]) and max(abs(psiIm[{1, 2, 4, 5}])) < 1.0e-5 * abs(psiRe[3]),
      "Schwarzschild Psi0, Psi1, Psi3, Psi4 are not zero in the radial tetrad");
    (IRe, IIm, JRe, JIm) := GeneralRelativity.Frames.weylInvariants(psiRe, psiIm);
    assert(abs(16.0 * IRe / inv.weylSquared - 1.0) < 1.0e-5, "16 Re I is not C_abcd C^abcd for Schwarzschild");
    assert(abs(JRe / (-psiRe[3]^3) - 1.0) < 1.0e-5, "Schwarzschild J is not -Psi2^3");
    assert(GeneralRelativity.Frames.classifyPetrovCurvature(bh) == GeneralRelativity.Frames.PetrovType.D, "Schwarzschild is not Petrov type D");
    cotton := GeneralRelativity.Metrics.cottonAt(GeneralRelativity.Metrics.schwarzschild, x, {rs}, step);
    // Cotton noise is round-off of nested differences, ~eps g / step^3; compare with the size of nabla Riemann ~ |R| / r.
    assert(maxAbs3(cotton) < 1.0e-5 * maxAbs4(bh.riemann) / r, "Schwarzschild Cotton tensor is not zero");

    // --- Kerr: type D, I^3 = 27 J^2, Weyl^2 = 16 Re I. ---
    kerr := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.kerr, x, {rs, 0.4 * rs}, step);
    C := GeneralRelativity.Tensors.weyl(kerr.g, kerr.riemann, kerr.ricci, kerr.ricciScalar);
    e := GeneralRelativity.Frames.orthonormalTetrad(kerr.g);
    assert(maxAbs4(GeneralRelativity.Frames.frameComponents4(C - kerr.riemann, e)) < 1.0e-5 * maxAbs4(GeneralRelativity.Frames.frameComponents4(kerr.riemann, e)),
      "Kerr Weyl tensor differs from Riemann");
    assert(tetradError(kerr.g, e) < 1.0e-12, "Kerr ZAMO tetrad is not orthonormal");
    (l, n, mRe, mIm) := GeneralRelativity.Frames.nullTetrad(e);
    (psiRe, psiIm) := GeneralRelativity.Frames.weylScalars(C, l, n, mRe, mIm);
    (IRe, IIm, JRe, JIm) := GeneralRelativity.Frames.weylInvariants(psiRe, psiIm);
    I3Re := IRe^3 - 3.0 * IRe * IIm^2;
    I3Im := 3.0 * IRe^2 * IIm - IIm^3;
    J2Re := JRe^2 - JIm^2;
    J2Im := 2.0 * JRe * JIm;
    assert(sqrt((I3Re - 27.0 * J2Re)^2 + (I3Im - 27.0 * J2Im)^2) < 1.0e-4 * sqrt(I3Re^2 + I3Im^2), "Kerr violates I^3 = 27 J^2");
    assert(abs(IIm) > 1.0e-3 * abs(IRe), "Kerr I should be complex (rotation makes Psi2 complex)");
    inv := GeneralRelativity.Tensors.curvatureInvariants(kerr);
    assert(abs(16.0 * IRe / inv.weylSquared - 1.0) < 1.0e-5, "16 Re I is not C_abcd C^abcd for Kerr");
    assert(GeneralRelativity.Frames.classifyPetrovCurvature(kerr) == GeneralRelativity.Frames.PetrovType.D, "Kerr is not Petrov type D");
    (E3, B3) := GeneralRelativity.Frames.weylElectricMagnetic(kerr.g, C);
    assert(abs(E3[1, 1] + E3[2, 2] + E3[3, 3]) < 1.0e-6 * Tests.Assertions.maxAbsMatrix(E3)
      and abs(B3[1, 1] + B3[2, 2] + B3[3, 3]) < 1.0e-6 * Tests.Assertions.maxAbsMatrix(E3), "Weyl E and B are not traceless");
    assert(Tests.Assertions.maxAbsMatrix(B3 - transpose(B3)) < 1.0e-6 * Tests.Assertions.maxAbsMatrix(E3), "Weyl magnetic part is not symmetric");
    assert(Tests.Assertions.maxAbsMatrix(B3) > 1.0e-3 * Tests.Assertions.maxAbsMatrix(E3), "Kerr has no gravitomagnetic Weyl part");

    // --- Tetrad inside the Kerr ergosphere (equator, r = 0.9 rs < r_ergo = rs, above r+ = 0.718 rs). ---
    g := GeneralRelativity.Metrics.kerr({0.0, 0.9 * rs, pi / 2.0, 0.0}, {rs, 0.45 * rs});
    (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(g);
    assert(g[1, 1] > 0.0 and gInv[1, 1] < 0.0, "Test point is not inside the ergosphere and outside the horizon");
    e := GeneralRelativity.Frames.orthonormalTetrad(g);
    assert(tetradError(g, e) < 1.0e-12, "ZAMO tetrad inside the ergosphere is not orthonormal");
    assert(e[1, 1] > 0.0 and e[1, 4] > 0.0, "ZAMO is not future-directed and co-rotating inside the ergosphere");
    // Observer tetrad of a boosted observer in Minkowski space.
    u := GeneralRelativity.Geodesics.fourVelocity(eta, {0.3, -0.4, 0.5}, 1.0);
    e := GeneralRelativity.Frames.observerTetrad(eta, u);
    assert(tetradError(eta, e) < 1.0e-12 and Tests.Assertions.maxAbsVector(e[1, :] - u) < 1.0e-12, "observerTetrad is not orthonormal with e0 = u");

    // --- Non-vacuum static sphere: quadratic identities, Weyl traceless, Petrov D, Cotton traces. ---
    x := {0.0, 1.0, 1.1, 0.4};
    star := GeneralRelativity.Metrics.curvatureAt(staticSphere, x, {0.0}, step);
    inv := GeneralRelativity.Tensors.curvatureInvariants(star);
    assert(abs(inv.ricciScalar) > 0.1 and inv.weylSquared > 1.0e-3, "Static test metric is unexpectedly vacuum or conformally flat");
    assert(abs(inv.weylSquared - (inv.kretschmann - 2.0 * inv.ricciSquared + inv.ricciScalar^2 / 3.0)) < 1.0e-8 * inv.kretschmann,
      "Weyl^2 = K - 2 Ric^2 + R^2/3 fails");
    assert(abs(inv.gaussBonnet - (inv.weylSquared - 2.0 * inv.ricciSquared + 2.0 * inv.ricciScalar^2 / 3.0)) < 1.0e-8 * inv.kretschmann,
      "Gauss-Bonnet = Weyl^2 - 2 Ric^2 + 2 R^2/3 fails");
    C := GeneralRelativity.Tensors.weyl(star.g, star.riemann, star.ricci, star.ricciScalar);
    assert(weylTraceError(star.gInv, C) < 1.0e-12 * maxAbs4(C) + 1.0e-12, "Non-vacuum Weyl tensor is not traceless");
    assert(abs(inv.cubic[3] - sum((star.gInv[i, i] * star.ricci[i, i])^3 for i in 1:4)) < 1.0e-10 * abs(inv.cubic[3]),
      "R_a^b R_b^c R_c^a disagrees with the diagonal Ricci eigenvalues");
    assert(abs(inv.cubic[1] - inv.ricciScalar^3) < 1.0e-12 * abs(inv.cubic[1]) and abs(inv.cubic[2] - inv.ricciScalar * inv.ricciSquared) < 1.0e-10 * abs(inv.cubic[2])
      and abs(inv.cubic[5] - inv.ricciScalar * inv.kretschmann) < 1.0e-10 * abs(inv.cubic[5]), "Products in the cubic basis are inconsistent");
    assert(GeneralRelativity.Frames.classifyPetrovCurvature(star) == GeneralRelativity.Frames.PetrovType.D, "Static spherical metric is not Petrov type D");
    cotton := GeneralRelativity.Metrics.cottonAt(staticSphere, x, {0.0}, step);
    assert(maxAbs3(cotton) > 1.0e-2, "Static test metric has an unexpectedly small Cotton tensor");
    cottonTrace := 0.0;
    for c in 1:4 loop
      cottonTrace := max(cottonTrace, abs(sum(star.gInv[a, b] * cotton[a, b, c] for a in 1:4, b in 1:4)));
      cottonTrace := max(cottonTrace, abs(sum(star.gInv[a, b] * cotton[a, c, b] for a in 1:4, b in 1:4)));
    end for;
    assert(cottonTrace < 1.0e-5 * maxAbs3(cotton), "Cotton tensor traces do not vanish (contracted Bianchi identity)");
    assert(max(abs(cotton[a, b, c] + cotton[a, c, b]) for a in 1:4, b in 1:4, c in 1:4) < 1.0e-14, "Cotton tensor is not antisymmetric in its last pair");

    // --- Schwarzschild-de Sitter: Einstein space, closed-form quadratic and cubic invariants. ---
    r := 6.0;
    x := {0.0, r, 1.1, 0.4};
    sds := GeneralRelativity.Metrics.curvatureAt(schwarzschildDeSitter, x, {rs, Lambda}, step);
    inv := GeneralRelativity.Tensors.curvatureInvariants(sds);
    m := rs / (2.0 * r^3);
    e := GeneralRelativity.Frames.orthonormalTetrad(sds.g);
    assert(Tests.Assertions.maxAbsMatrix(GeneralRelativity.Frames.frameComponents2(GeneralRelativity.Tensors.traceFreeRicci(sds.g, sds.ricci, sds.ricciScalar), e)) < 1.0e-5 * Lambda,
      "Schwarzschild-de Sitter is not an Einstein space");
    assert(abs(inv.ricciScalar / (4.0 * Lambda) - 1.0) < 1.0e-5, "Schwarzschild-de Sitter R is not 4 Lambda");
    assert(abs(inv.weylSquared / (48.0 * m^2) - 1.0) < 1.0e-5, "Schwarzschild-de Sitter Weyl^2 is not 48 M^2/r^6");
    assert(abs(inv.kretschmann / (48.0 * m^2 + 8.0 * Lambda^2 / 3.0) - 1.0) < 1.0e-5, "Schwarzschild-de Sitter K is not 48 M^2/r^6 + 8 Lambda^2/3");
    assert(abs(inv.gaussBonnet / (48.0 * m^2 + 8.0 * Lambda^2 / 3.0) - 1.0) < 1.0e-5, "Einstein-space Gauss-Bonnet density is not K");
    assert(abs(inv.cubic[1] / (64.0 * Lambda^3) - 1.0) < 1.0e-4 and abs(inv.cubic[2] / (16.0 * Lambda^3) - 1.0) < 1.0e-4
      and abs(inv.cubic[3] / (4.0 * Lambda^3) - 1.0) < 1.0e-4 and abs(inv.cubic[4] / (4.0 * Lambda^3) - 1.0) < 1.0e-4,
      "Einstein-space Ricci cubic invariants are not 64, 16, 4, 4 Lambda^3");
    // R_ab^cd = C_ab^cd + (Lambda/3)(delta delta - delta delta): tr M^3 = 8 sum (c_i + Lambda/3)^3 over the six bivectors.
    expected := 96.0 * m^3 + 96.0 * Lambda * m^2 + 16.0 * Lambda^3 / 9.0;
    assert(abs(inv.cubic[6] / expected - 1.0) < 1.0e-5, "Schwarzschild-de Sitter Riemann cubed disagrees with the closed form");
    assert(abs(inv.weylCubed / (96.0 * m^3) - 1.0) < 1.0e-5, "Schwarzschild-de Sitter Weyl cubed is not 96 M^3/r^9");
    cotton := GeneralRelativity.Metrics.cottonAt(schwarzschildDeSitter, x, {rs, Lambda}, step);
    assert(maxAbs3(cotton) < 1.0e-5 * maxAbs4(sds.riemann) / r, "Einstein-space Cotton tensor is not zero");

    // --- Conformally flat metric: Weyl and Cotton vanish, Schouten does not, Petrov O. ---
    x := {0.3, 0.5, -0.4, 0.8};
    conf := GeneralRelativity.Metrics.curvatureAt(conformallyFlat, x, {0.0}, step);
    C := GeneralRelativity.Tensors.weyl(conf.g, conf.riemann, conf.ricci, conf.ricciScalar);
    P := GeneralRelativity.Tensors.schouten(conf.g, conf.ricci, conf.ricciScalar);
    assert(Tests.Assertions.maxAbsMatrix(P) > 1.0e-2, "Conformally flat test metric has a negligible Schouten tensor");
    assert(maxAbs4(C) < 1.0e-6 * maxAbs4(conf.riemann), "Weyl tensor of a conformally flat metric is not zero");
    assert(GeneralRelativity.Frames.classifyPetrovCurvature(conf) == GeneralRelativity.Frames.PetrovType.O, "Conformally flat metric is not Petrov type O");
    cotton := GeneralRelativity.Metrics.cottonAt(conformallyFlat, x, {0.0}, step);
    assert(maxAbs3(cotton) < 1.0e-5 * maxAbs4(conf.riemann), "Cotton tensor of a conformally flat metric is not zero");

    // --- pp-wave: vacuum, type N, both Weyl invariants zero (in any tetrad). ---
    x := {0.2, 0.7, 0.3, -0.1};
    wave := GeneralRelativity.Metrics.curvatureAt(ppWave, x, {0.0}, step);
    assert(Tests.Assertions.maxAbsMatrix(wave.ricci) < 1.0e-6, "pp-wave is not vacuum");
    assert(maxAbs4(wave.riemann) > 0.1, "pp-wave curvature is unexpectedly small");
    assert(abs(wave.kretschmann) < 1.0e-6, "pp-wave Kretschmann scalar is not zero");
    assert(GeneralRelativity.Frames.classifyPetrovCurvature(wave) == GeneralRelativity.Frames.PetrovType.N, "pp-wave is not Petrov type N");
    C := GeneralRelativity.Tensors.weyl(wave.g, wave.riemann, wave.ricci, wave.ricciScalar);
    e := GeneralRelativity.Frames.orthonormalTetrad(wave.g);
    (l, n, mRe, mIm) := GeneralRelativity.Frames.nullTetrad(e);
    (psiRe, psiIm) := GeneralRelativity.Frames.weylScalars(C, l, n, mRe, mIm);
    (IRe, IIm, JRe, JIm) := GeneralRelativity.Frames.weylInvariants(psiRe, psiIm);
    lambda := max(max(abs(psiRe)), max(abs(psiIm)));
    assert(lambda > 0.1 and abs(IRe) + abs(IIm) < 1.0e-5 * lambda^2 and abs(JRe) + abs(JIm) < 1.0e-5 * lambda^3,
      "pp-wave invariants I and J are not zero");

    // --- Generic metric: type I, and I^3 != 27 J^2. ---
    x := {0.1, 0.2, -0.1, 0.15};
    generic := GeneralRelativity.Metrics.curvatureAt(genericMetric, x, {0.0}, step);
    assert(GeneralRelativity.Frames.classifyPetrovCurvature(generic, 1.0e-2) == GeneralRelativity.Frames.PetrovType.I, "Generic metric is not Petrov type I");
    C := GeneralRelativity.Tensors.weyl(generic.g, generic.riemann, generic.ricci, generic.ricciScalar);
    assert(weylTraceError(generic.gInv, C) < 1.0e-10 * maxAbs4(C), "Generic Weyl tensor is not traceless");
    inv := GeneralRelativity.Tensors.curvatureInvariants(generic);
    assert(abs(inv.weylSquared - (inv.kretschmann - 2.0 * inv.ricciSquared + inv.ricciScalar^2 / 3.0)) < 1.0e-8 * maxAbs4(generic.riemann)^2,
      "Weyl^2 identity fails for the generic metric");

    // --- Petrov decision tree on hand-made Q = E + iB matrices (x = (1, i, 0) is a complex null vector). ---
    lambda := 0.3;
    assert(GeneralRelativity.Frames.classifyPetrovMatrix(diagonal({1.0, 2.0, -3.0}), zeros(3, 3), 1.0) == GeneralRelativity.Frames.PetrovType.I, "Matrix type I");
    assert(GeneralRelativity.Frames.classifyPetrovMatrix(lambda * diagonal({1.0, 1.0, -2.0}), zeros(3, 3), 1.0) == GeneralRelativity.Frames.PetrovType.D, "Matrix type D");
    assert(GeneralRelativity.Frames.classifyPetrovMatrix(lambda * diagonal({1.0, 1.0, -2.0}) + [1, 0, 0; 0, -1, 0; 0, 0, 0], [0, 1, 0; 1, 0, 0; 0, 0, 0], 1.0)
      == GeneralRelativity.Frames.PetrovType.II, "Matrix type II");
    assert(GeneralRelativity.Frames.classifyPetrovMatrix([0, 0, 1; 0, 0, 0; 1, 0, 0], [0, 0, 0; 0, 0, 1; 0, 1, 0], 1.0) == GeneralRelativity.Frames.PetrovType.III, "Matrix type III");
    assert(GeneralRelativity.Frames.classifyPetrovMatrix([1, 0, 0; 0, -1, 0; 0, 0, 0], [0, 1, 0; 1, 0, 0; 0, 0, 0], 1.0) == GeneralRelativity.Frames.PetrovType.N, "Matrix type N");
    assert(GeneralRelativity.Frames.classifyPetrovMatrix(1.0e-9 * diagonal({1.0, 2.0, -3.0}), zeros(3, 3), 1.0) == GeneralRelativity.Frames.PetrovType.O, "Matrix type O");

    // --- Energy conditions. ---
    // Perfect fluid at rest.
    T := GeneralRelativity.StressEnergy.perfectFluid(eta, {1.0, 0.0, 0.0, 0.0}, 3.0, 0.5, 1.0);
    (nec, wec, sec, dec, rho, pr, typeI, u) := GeneralRelativity.Frames.energyConditions(eta, T);
    assert(typeI and abs(rho - 3.0) < 1.0e-12 and Tests.Assertions.maxAbsVector(pr - fill(0.5, 3)) < 1.0e-12, "Fluid at rest: rho, p wrong");
    assert(nec and wec and sec and dec, "Ordinary fluid violates an energy condition");
    // Boosted perfect fluid: same rho and p, rest frame recovered.
    u := GeneralRelativity.Geodesics.fourVelocity(eta, {0.6, -0.3, 0.2}, 1.0);
    T := GeneralRelativity.StressEnergy.perfectFluid(eta, u, 3.0, 0.5, 1.0);
    (nec, wec, sec, dec, rho, pr, typeI, k) := GeneralRelativity.Frames.energyConditions(eta, T);
    assert(typeI and abs(rho - 3.0) < 1.0e-10 and Tests.Assertions.maxAbsVector(pr - fill(0.5, 3)) < 1.0e-10, "Boosted fluid: rho, p wrong");
    assert(Tests.Assertions.maxAbsVector(k - u) < 1.0e-8, "Boosted fluid: rest-frame four-velocity not recovered");
    assert(nec and wec and sec and dec, "Boosted ordinary fluid violates an energy condition");
    // Fluid moving through the Kerr metric (non-diagonal g, non-ZAMO u).
    g := kerr.g;
    u := {1.0, 0.05, -0.01, 0.01};
    u := u / sqrt(-(u * g * u));
    T := GeneralRelativity.StressEnergy.perfectFluid(g, u, 2.0, 0.3, 1.0);
    (nec, wec, sec, dec, rho, pr, typeI, k) := GeneralRelativity.Frames.energyConditions(g, T);
    assert(typeI and abs(rho - 2.0) < 1.0e-10 and Tests.Assertions.maxAbsVector(pr - fill(0.3, 3)) < 1.0e-10
      and Tests.Assertions.maxAbsVector(k - u) < 1.0e-8, "Fluid in Kerr: rho, p or u wrong");
    // Anisotropic stress in a rotated frame: principal pressures recovered.
    T := diagonal({1.0, -0.2, 0.4, 0.9});
    T[3, 4] := 0.25;
    T[4, 3] := 0.25;
    (nec, wec, sec, dec, rho, pr, typeI, k) := GeneralRelativity.Frames.energyConditions(eta, T);
    assert(typeI and abs(rho - 1.0) < 1.0e-12 and abs(pr[1] + 0.2) < 1.0e-12
      and abs(pr[2] - (0.65 - sqrt(0.0625 + 0.0625))) < 1.0e-12 and abs(pr[3] - (0.65 + sqrt(0.125))) < 1.0e-12,
      "Anisotropic stress: principal pressures wrong");
    // Radiation.
    T := GeneralRelativity.StressEnergy.radiation(eta, GeneralRelativity.Geodesics.fourVelocity(eta, {0.0, 0.5, 0.0}, 1.0), 3.0, 1.0);
    (nec, wec, sec, dec, rho, pr, typeI, k) := GeneralRelativity.Frames.energyConditions(eta, T);
    assert(typeI and abs(rho - 3.0) < 1.0e-10 and Tests.Assertions.maxAbsVector(pr - fill(1.0, 3)) < 1.0e-10, "Radiation: rho, p wrong");
    assert(nec and wec and sec and dec, "Radiation violates an energy condition");
    // Positive cosmological constant: rho = -p > 0; NEC marginal, WEC and DEC hold, SEC fails.
    T := GeneralRelativity.StressEnergy.cosmologicalConstant(kerr.g, Lambda, 1.0, 1.0);
    (nec, wec, sec, dec, rho, pr, typeI, k) := GeneralRelativity.Frames.energyConditions(kerr.g, T);
    assert(typeI and abs(rho / (Lambda / (8.0 * pi)) - 1.0) < 1.0e-10 and Tests.Assertions.maxAbsVector(pr + fill(rho, 3)) < 1.0e-10 * rho,
      "Cosmological constant: rho, p wrong");
    assert(nec and wec and dec and not sec, "Cosmological constant: expected NEC, WEC, DEC true and SEC false");
    // Negative energy density with positive pressure: NEC and SEC hold, WEC and DEC fail.
    T := GeneralRelativity.StressEnergy.perfectFluid(eta, GeneralRelativity.Geodesics.fourVelocity(eta, {0.0, 0.0, 0.4}, 1.0), -0.2, 0.5, 1.0);
    (nec, wec, sec, dec, rho, pr, typeI, k) := GeneralRelativity.Frames.energyConditions(eta, T);
    assert(typeI and abs(rho + 0.2) < 1.0e-10 and Tests.Assertions.maxAbsVector(pr - fill(0.5, 3)) < 1.0e-10, "Negative-energy fluid: rho, p wrong");
    assert(nec and not wec and sec and not dec, "Negative-energy fluid: expected NEC, SEC true and WEC, DEC false");
    // Phantom fluid rho + p < 0: every condition fails.
    T := GeneralRelativity.StressEnergy.perfectFluid(eta, {1.0, 0.0, 0.0, 0.0}, 1.0, -1.5, 1.0);
    (nec, wec, sec, dec, rho, pr, typeI, k) := GeneralRelativity.Frames.energyConditions(eta, T);
    assert(typeI and not (nec or wec or sec or dec), "Phantom fluid: every energy condition should fail");
    // Null dust T_ab = k_a k_b is Hawking-Ellis type II, not type I.
    k := {-1.0, 1.0, 0.0, 0.0};
    T := {k[a] * k[b] for b in 1:4, a in 1:4};
    (nec, wec, sec, dec, rho, pr, typeI, u) := GeneralRelativity.Frames.energyConditions(eta, T);
    assert(not typeI, "Null dust was reported as type I");

    passed := true;
  end run;

  parameter Boolean passed = run();
equation
  assert(passed, "GeneralRelativity curvature assertions did not complete");
  annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 1.0e-8, Interval = 0.1));
end GeneralRelativityCurvatureTests;
