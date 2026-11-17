within Tests;
model GeneralRelativityIndexTests "Inline index-notation evaluator: contractions, raising and lowering, sums, fixed components, user tensors and syntax variants"
  function staticSphere "Non-vacuum static spherically symmetric test metric in {ct, r, theta, phi}"
    extends GeneralRelativity.Interfaces.MetricFunction;
  algorithm
    g := zeros(4, 4);
    g[1, 1] := -(1.0 + x[2]^2);
    g[2, 2] := 1.0 / (1.0 - x[2]^2 / 4.0);
    g[3, 3] := x[2]^2;
    g[4, 4] := x[2]^2 * sin(x[3])^2;
  end staticSphere;

  function maxAbs3 "Maximum absolute element of a rank-3 array"
    input Real A[:, :, :];
    output Real m;
  algorithm
    m := max(abs(A[a, b, c]) for a in 1:size(A, 1), b in 1:size(A, 2), c in 1:size(A, 3));
  end maxAbs3;

  function maxAbs4 "Maximum absolute element of a rank-4 array"
    input Real A[:, :, :, :];
    output Real m;
  algorithm
    m := max(abs(A[a, b, c, d]) for a in 1:size(A, 1), b in 1:size(A, 2), c in 1:size(A, 3), d in 1:size(A, 4));
  end maxAbs4;

  function run
    input Real t "Unused; a time argument keeps the call out of the front end's constant evaluation";
    output Boolean passed;
  protected
    import GeneralRelativity.IndexNotation.scalar;
    import GeneralRelativity.IndexNotation.vector;
    import GeneralRelativity.IndexNotation.matrix;
    import GeneralRelativity.IndexNotation.tensor3;
    import GeneralRelativity.IndexNotation.tensor4;
    import GeneralRelativity.IndexNotation.evaluate;
    import Tests.Assertions.maxAbsMatrix;
    import Tests.Assertions.maxAbsVector;
    constant Real pi = 3.141592653589793;
    constant Real step[4] = {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3};
    constant Real eta[4, 4] = [-1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
    GeneralRelativity.Tensors.Curvature bh;
    GeneralRelativity.Tensors.Curvature star;
    GeneralRelativity.Tensors.Invariants inv;
    GeneralRelativity.IndexNotation.TensorSystem sys;
    GeneralRelativity.IndexNotation.TensorValue v;
    Real riemannScale;
    Real ricciScale;
    Real u[4];
    Real w[4];
    Real T[4, 4];
    Real det;
  algorithm
    // ------------------------------------------------------------------
    // Schwarzschild, vacuum: metric identities and the Kretschmann scalar.
    bh := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.schwarzschild, {0.0, 10.0, pi / 2.0, 0.0}, {2.0}, step);
    sys := GeneralRelativity.IndexNotation.systemFromCurvature(bh);
    riemannScale := maxAbs4(bh.riemann);
    assert(abs(scalar(sys, "g_{mu nu} g^{mu nu}") - 4.0) < 1.0e-12, "g_{mu nu} g^{mu nu} is not 4");
    assert(abs(scalar(sys, "delta^a_a") - 4.0) < 1.0e-12, "delta^a_a is not 4");
    assert(abs(scalar(sys, "g_ab g^ab") - 4.0) < 1.0e-12, "bare index groups g_ab g^ab failed");
    assert(abs(scalar(sys, "g_{\\mu\\nu} g^{\\mu \\nu}") - 4.0) < 1.0e-12, "LaTeX Greek indices failed");
    assert(abs(scalar(sys, "g_{mu,nu}*g^{mu,nu}") - 4.0) < 1.0e-12, "comma-separated indices with * failed");
    assert(abs(scalar(sys, "  g_{ab}g^{ab}  ") - 4.0) < 1.0e-12, "whitespace variant failed");
    assert(maxAbsMatrix(matrix(sys, "g^{ab} g_{bc}") - identity(4)) < 1.0e-12, "g^{ab} g_{bc} is not the identity");
    assert(maxAbsMatrix(matrix(sys, "g^{ab}") - bh.gInv) < 1.0e-12, "g^{ab} is not the inverse metric");
    assert(maxAbsMatrix(matrix(sys, "g_{ab}") - bh.g) < 1.0e-12, "g_{ab} is not the metric");
    assert(abs(scalar(sys, "R_{abcd} R^{abcd}") / bh.kretschmann - 1.0) < 1.0e-10, "R_{abcd} R^{abcd} is not the Kretschmann scalar");
    assert(abs(scalar(sys, "R_{abcd} R_{abcd}") / bh.kretschmann - 1.0) < 1.0e-10, "same-position contraction did not insert the metric");
    assert(abs(scalar(sys, "R^{abcd} R^{abcd}") / bh.kretschmann - 1.0) < 1.0e-10, "both-upper contraction did not insert the metric");
    assert(abs(scalar(sys, "K") / bh.kretschmann - 1.0) < 1.0e-12, "K is not the Kretschmann scalar");
    assert(abs(scalar(sys, "R_{abcd} R^{abcd} / K") - 1.0) < 1.0e-10, "division by a scalar factor failed");
    assert(abs(scalar(sys, "R_{abcd} R^{abcd} / (12 * 4 / 1e6)") - 1.0) < 1.0e-4, "closed-form Kretschmann 12 rs^2 / r^6 failed");
    assert(abs(scalar(sys, "R")) < 1.0e-6 * riemannScale, "Schwarzschild Ricci scalar is not zero");
    assert(maxAbsMatrix(matrix(sys, "R^{a}_{bac}")) < 1.0e-6 * riemannScale, "Schwarzschild Ricci tensor from the Riemann contraction is not zero");
    assert(maxAbs4(tensor4(sys, "C_{abcd}") - GeneralRelativity.Tensors.weyl(bh.g, bh.riemann, bh.ricci, bh.ricciScalar)) < 1.0e-16, "C_{abcd} is not the Weyl tensor");
    // Finite-difference Ricci noise (about 1e-7) enters C = R - g wedge P multiplied by metric components of order r^2.
    assert(maxAbs4(tensor4(sys, "C_{abcd} - R_{abcd}")) < 1.0e-4, "Weyl tensor differs from Riemann in vacuum");
    assert(maxAbs4(tensor4(sys, "R_{abcd}") - bh.riemann) < 1.0e-14, "R_{abcd} does not reproduce the Riemann array");
    // Raising then lowering again returns the array; free index order follows first appearance: c d e f.
    assert(maxAbs4(tensor4(sys, "R^{ab}_{cd} g_{ae} g_{bf}") - bh.riemann) < 1.0e-12 * riemannScale, "R^{ab}_{cd} g_{ae} g_{bf} is not R_{cdef} = R_{efcd}");
    // Symmetries of the Riemann tensor.
    assert(maxAbs4(tensor4(sys, "R_{abcd} + R_{bacd}")) < 1.0e-12 * riemannScale, "R_{abcd} is not antisymmetric in ab");
    assert(maxAbs4(tensor4(sys, "R_{abcd} - R_{cdab}")) < 1.0e-12 * riemannScale, "R_{abcd} lacks pair symmetry");
    assert(maxAbs4(tensor4(sys, "R_{abcd} + R_{acdb} + R_{adbc}")) < 1.0e-12 * riemannScale, "first Bianchi identity failed");
    // Christoffel symbols: stored lowered, so raising the first index gives back the array.
    assert(maxAbs3(tensor3(sys, "Gamma^i_{jk}") - bh.Gamma) < 1.0e-12, "Gamma^i_{jk} does not reproduce the Christoffel array");
    assert(maxAbs3(tensor3(sys, "Gamma^{i}_{jk} - Gamma^{i}_{kj}")) < 1.0e-14, "Christoffel symbols are not symmetric in the lower pair");
    w := {sum(bh.Gamma[i, i, j] for i in 1:4) for j in 1:4};
    assert(maxAbsVector(vector(sys, "Gamma^i_{ij}") - w) < 1.0e-12, "Gamma^i_{ij} contraction failed");
    // Levi-Civita tensor.
    det := Modelica.Math.Matrices.det(bh.g);
    assert(abs(scalar(sys, "epsilon_{abcd} epsilon^{abcd}") + 24.0) < 1.0e-10, "epsilon_{abcd} epsilon^{abcd} is not -24");
    assert(abs(scalar(sys, "eps_{0123}") - sqrt(-det)) < 1.0e-12 * sqrt(-det), "epsilon_{0123} is not sqrt(-det g)");
    assert(abs(scalar(sys, "epsilon_{1023}") + sqrt(-det)) < 1.0e-12 * sqrt(-det), "epsilon_{1023} is not -sqrt(-det g)");
    // Fixed components in physics labels 0..3, before and after raising.
    assert(abs(scalar(sys, "g_{00}") - bh.g[1, 1]) < 1.0e-14, "g_{00} is not g[1, 1]");
    assert(abs(scalar(sys, "g^{11}") - bh.gInv[2, 2]) < 1.0e-14, "g^{11} is not gInv[2, 2]");
    assert(abs(scalar(sys, "R_{0101}") - bh.riemann[1, 2, 1, 2]) < 1.0e-16, "R_{0101} is not riemann[1, 2, 1, 2]");
    assert(abs(scalar(sys, "R^{0}_{101}") - bh.gInv[1, 1] * bh.riemann[1, 2, 1, 2]) < 1.0e-16, "R^{0}_{101} is not gInv[1, 1] riemann[1, 2, 1, 2]");
    assert(abs(scalar(sys, "g_{0 mu} g^{mu 0}") - 1.0) < 1.0e-12, "g_{0 mu} g^{mu 0} is not delta^0_0 = 1");
    assert(abs(scalar(sys, "R_{0 a 0 b} g^{ab}") - sum(bh.gInv[a, b] * bh.riemann[1, a, 1, b] for a in 1:4, b in 1:4)) < 1.0e-16, "fixed indices with a full contraction failed");
    assert(maxAbsVector(vector(sys, "R_{0a0b} g^{1b}") - {sum(bh.gInv[2, b] * bh.riemann[1, a, 1, b] for b in 1:4) for a in 1:4}) < 1.0e-16, "mixed fixed and free indices with contraction failed");
    assert(maxAbsVector(vector(sys, "g_{1a}") - bh.g[2, :]) < 1.0e-14, "g_{1a} is not the second row of g");
    // ------------------------------------------------------------------
    // Non-vacuum static metric: Ricci, Einstein, Weyl, Schouten, trace-free Ricci, stress-energy.
    star := GeneralRelativity.Metrics.curvatureAt(staticSphere, {0.0, 0.7, 1.1, 0.3}, {0.0}, step);
    inv := GeneralRelativity.Tensors.curvatureInvariants(star);
    sys := GeneralRelativity.IndexNotation.systemFromCurvature(star, 1.0, 1.0, 0.0);
    riemannScale := maxAbs4(star.riemann);
    ricciScale := maxAbsMatrix(star.ricci);
    assert(ricciScale > 1.0e-3, "Test metric is not curved enough to be a useful check");
    assert(maxAbsMatrix(matrix(sys, "R^{a}_{bac}") - star.ricci) < 1.0e-10 * ricciScale, "R^{a}_{bac} is not the Ricci tensor");
    assert(maxAbsMatrix(matrix(sys, "g^{ac} R_{abcd}") - star.ricci) < 1.0e-10 * ricciScale, "g^{ac} R_{abcd} is not the Ricci tensor");
    assert(maxAbsMatrix(matrix(sys, "Ric_{ab}") - star.ricci) < 1.0e-14, "Ric alias failed");
    assert(maxAbsMatrix(matrix(sys, "Ricci_{ab}") - matrix(sys, "R_{ab}")) < 1.0e-16, "Ricci alias failed");
    assert(abs(scalar(sys, "g^{ab} R_{ab}") / star.ricciScalar - 1.0) < 1.0e-10, "g^{ab} R_{ab} is not the Ricci scalar");
    assert(abs(scalar(sys, "R^{a}_{a}") / star.ricciScalar - 1.0) < 1.0e-10, "R^a_a is not the Ricci scalar");
    assert(abs(scalar(sys, "R") / star.ricciScalar - 1.0) < 1.0e-14, "R with no indices is not the Ricci scalar");
    assert(maxAbsMatrix(matrix(sys, "R_{ab} - 1/2 g_{ab} R - G_{ab}")) < 1.0e-12 * ricciScale, "Einstein identity failed");
    assert(maxAbsMatrix(matrix(sys, "Einstein_{ab} - G_{ab}")) < 1.0e-16, "Einstein alias failed");
    assert(abs(scalar(sys, "G^{a}_{a}") / star.ricciScalar + 1.0) < 1.0e-10, "G^a_a is not -R");
    assert(abs(scalar(sys, "(R_{ab} - G_{ab}) g^{ab}") / star.ricciScalar - 2.0) < 1.0e-10, "(R_{ab} - G_{ab}) g^{ab} is not 2R");
    assert(abs(scalar(sys, "R_{ab} R^{ab}") / inv.ricciSquared - 1.0) < 1.0e-10, "R_{ab} R^{ab} is not the Ricci-squared invariant");
    assert(abs(scalar(sys, "C_{abcd} C^{abcd}") / inv.weylSquared - 1.0) < 1.0e-8, "C_{abcd} C^{abcd} is not the Weyl invariant");
    assert(maxAbsMatrix(matrix(sys, "C^{a}_{bad}")) < 1.0e-10 * riemannScale, "Weyl tensor is not trace free");
    assert(maxAbs4(tensor4(sys, "Weyl_{abcd} - C_{abcd}")) < 1.0e-16, "Weyl alias failed");
    assert(maxAbsMatrix(matrix(sys, "P_{ab}") - GeneralRelativity.Tensors.schouten(star.g, star.ricci, star.ricciScalar)) < 1.0e-16, "Schouten tensor missing");
    assert(maxAbsMatrix(matrix(sys, "Schouten_{ab} - P_{ab}")) < 1.0e-16, "Schouten alias failed");
    assert(maxAbsMatrix(matrix(sys, "S_{ab} - R_{ab} + 1/4 g_{ab} R")) < 1.0e-12 * ricciScale, "trace-free Ricci tensor is wrong");
    assert(abs(scalar(sys, "S^{a}_{a}")) < 1.0e-12 * ricciScale, "trace-free Ricci tensor has a trace");
    T := GeneralRelativity.Tensors.stressEnergyFromEinstein(star.einstein, star.g, 1.0, 1.0, 0.0);
    assert(maxAbsMatrix(matrix(sys, "T_{ab}") - T) < 1.0e-16, "stress-energy tensor missing");
    // Trace of the field equations in geometric units: 8 pi T = -R.
    assert(abs(scalar(sys, "8 * 3.141592653589793 T_{ab} g^{ab} + R")) < 1.0e-10 * abs(star.ricciScalar), "trace of the field equations failed");
    // Weyl decomposition: C = R - g wedge P (Kulkarni-Nomizu), written out in index notation.
    assert(maxAbs4(tensor4(sys, "C_{abcd} - R_{abcd} + g_{ac} P_{bd} - g_{ad} P_{bc} - g_{bc} P_{ad} + g_{bd} P_{ac}")) < 1.0e-10 * riemannScale,
      "Weyl decomposition R = C + g wedge P failed");
    // Sums with permuted labels and mixed positions.
    assert(maxAbsMatrix(matrix(sys, "R_{ab} - R_{ba}")) < 1.0e-14 * ricciScale, "R_{ab} - R_{ba} is not zero");
    assert(maxAbsMatrix(matrix(sys, "R^{a}_{b} - R_{b}^{a}")) < 1.0e-12 * ricciScale, "R^a_b - R_b^a is not zero");
    assert(maxAbsMatrix(matrix(sys, "R_{ab} - R^{ab}")) < 1.0e-12 * ricciScale, "adding terms in different positions did not lower the second");
    assert(maxAbsMatrix(matrix(sys, "0.5*R_{ab} + R_{ab}/2 - R_{ab}")) < 1.0e-14 * ricciScale, "scalar factors in a sum failed");
    assert(abs(scalar(sys, "-(g_{ab} g^{ab}) + 8")) - 4.0 < 1.0e-12, "leading minus with parentheses failed");
    assert(abs(scalar(sys, "2 g_{ab} g^{ab} / 4") - 2.0) < 1.0e-12, "2 g_{ab} g^{ab} / 4 is not 2");
    assert(abs(scalar(sys, "1e-3 g_{ab} g^{ab}") - 4.0e-3) < 1.0e-15, "exponent literal failed");
    assert(abs(scalar(sys, "0.5*R - R/2")) < 1.0e-14 * abs(star.ricciScalar), "0.5 R - R/2 is not zero");
    assert(abs(scalar(sys, "+R - R")) < 1.0e-16, "leading plus failed");
    // evaluate: free index order and positions.
    v := evaluate(sys, "R^{a}_{bc d}");
    assert(v.rank == 4 and v.index[1] == "a" and v.index[2] == "b" and v.index[3] == "c" and v.index[4] == "d", "free index labels are wrong");
    assert(v.up[1] and not v.up[2] and not v.up[3] and not v.up[4], "free index positions are wrong");
    u := {1.0 / sqrt(-star.g[1, 1]), 0.0, 0.0, 0.0};
    sys := GeneralRelativity.IndexNotation.defineVector(sys, "u", u);
    v := evaluate(sys, "G_{ab} u^b");
    assert(v.rank == 1 and v.index[1] == "a" and not v.up[1], "contraction with a user vector left the wrong free index");
    assert(maxAbsVector(v.comp[1:4] - star.einstein * u) < 1.0e-14, "G_{ab} u^b components are wrong");
    // toString covers every rank branch.
    assert(Modelica.Utilities.Strings.length(GeneralRelativity.IndexNotation.toString(evaluate(sys, "R"))) > 0, "toString of a scalar is empty");
    assert(Modelica.Utilities.Strings.length(GeneralRelativity.IndexNotation.toString(evaluate(sys, "G_{ab} u^b"))) > 0, "toString of a vector is empty");
    assert(Modelica.Utilities.Strings.length(GeneralRelativity.IndexNotation.toString(evaluate(sys, "G_{ab}"))) > 0, "toString of a matrix is empty");
    assert(Modelica.Utilities.Strings.length(GeneralRelativity.IndexNotation.toString(evaluate(sys, "Gamma^i_{jk}"))) > 0, "toString of a rank-3 tensor is empty");
    passed := true;
  end run;

  function runUserTensors "User-defined tensors on the Schwarzschild system and on flat space"
    input Real t "Unused; a time argument keeps the call out of the front end's constant evaluation";
    output Boolean passed;
  protected
    import GeneralRelativity.IndexNotation.scalar;
    import GeneralRelativity.IndexNotation.vector;
    import GeneralRelativity.IndexNotation.matrix;
    import Tests.Assertions.maxAbsMatrix;
    import Tests.Assertions.maxAbsVector;
    constant Real pi = 3.141592653589793;
    constant Real step[4] = {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3};
    constant Real eta[4, 4] = [-1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
    GeneralRelativity.Tensors.Curvature bh;
    GeneralRelativity.IndexNotation.TensorSystem sys;
    GeneralRelativity.IndexNotation.TensorSystem flat;
    Real u[4];
    Real w[4];
    Real F[4, 4];
    Real s;
  algorithm
    bh := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.schwarzschild, {0.0, 10.0, pi / 2.0, 0.0}, {2.0}, step);
    sys := GeneralRelativity.IndexNotation.systemFromCurvature(bh);
    // Static observer: contravariant four-velocity normalized to g(u, u) = -1.
    u := {1.0 / sqrt(-bh.g[1, 1]), 0.0, 0.0, 0.0};
    sys := GeneralRelativity.IndexNotation.defineVector(sys, "u", u);
    assert(abs(scalar(sys, "g_{ab} u^a u^b") + 1.0) < 1.0e-12, "g_{ab} u^a u^b is not -1");
    assert(abs(scalar(sys, "u_a u^a") + 1.0) < 1.0e-12, "u_a u^a is not -1");
    assert(abs(scalar(sys, "u^a u_a") + 1.0) < 1.0e-12, "u^a u_a is not -1");
    assert(maxAbsVector(vector(sys, "u^a") - u) < 1.0e-14, "u^a does not return the components as defined");
    assert(maxAbsVector(vector(sys, "u_a") - bh.g * u) < 1.0e-14, "u_a is not g_ab u^b");
    assert(maxAbsVector(vector(sys, "g_{ab} u^b") - bh.g * u) < 1.0e-14, "g_{ab} u^b is not the lowered vector");
    // Covariant definition round trip.
    sys := GeneralRelativity.IndexNotation.defineVector(sys, "p", bh.g * u, false);
    assert(abs(scalar(sys, "p_a u^a") + 1.0) < 1.0e-12, "p_a u^a is not -1 for p = lowered u");
    assert(maxAbsVector(vector(sys, "p^a") - u) < 1.0e-12, "p^a is not u^a");
    // Energy density seen by the static observer: T_ab u^a u^b, zero in vacuum.
    assert(abs(scalar(sys, "T_{ab} u^a u^b")) < 1.0e-6 * Tests.Assertions.maxAbsMatrix(bh.riemann[1, :, 1, :]) + 1.0e-12, "Schwarzschild energy density is not zero");
    // Tidal tensor R_{a c b d} u^c u^d is symmetric and orthogonal to u.
    assert(maxAbsMatrix(matrix(sys, "R_{acbd} u^c u^d - R_{bcad} u^c u^d")) < 1.0e-16, "tidal tensor is not symmetric");
    assert(maxAbsVector(vector(sys, "R_{acbd} u^c u^d u^a")) < 1.0e-16, "tidal tensor is not orthogonal to u");
    // Scalars scale vectors; redefinition replaces components.
    sys := GeneralRelativity.IndexNotation.defineScalar(sys, "m", 3.0);
    assert(maxAbsVector(vector(sys, "m u^a") - 3.0 * u) < 1.0e-14, "m u^a is not 3 u^a");
    sys := GeneralRelativity.IndexNotation.defineVector(sys, "u", 2.0 * u);
    assert(abs(scalar(sys, "u_a u^a") + 4.0) < 1.0e-12, "redefining u did not replace its components");
    // Antisymmetric field tensor on flat space: F_{ab} F^{ab} = 2 (B^2 - E^2).
    flat := GeneralRelativity.IndexNotation.systemFromMetric(eta);
    assert(abs(scalar(flat, "g_{ab} g^{ab}") - 4.0) < 1.0e-14, "flat metric trace is not 4");
    F := zeros(4, 4);
    F[1, 2] := 0.3;
    F[2, 1] := -0.3;
    F[1, 3] := -0.1;
    F[3, 1] := 0.1;
    F[2, 3] := 0.7;
    F[3, 2] := -0.7;
    F[3, 4] := 0.2;
    F[4, 3] := -0.2;
    flat := GeneralRelativity.IndexNotation.defineMatrix(flat, "F", F);
    s := 2.0 * (0.7^2 + 0.2^2 - 0.3^2 - 0.1^2);
    assert(abs(scalar(flat, "F_{ab} F^{ab}") - s) < 1.0e-12, "F_{ab} F^{ab} is not 2 (B^2 - E^2)");
    assert(maxAbsMatrix(matrix(flat, "F_{ab} + F_{ba}")) < 1.0e-16, "F is not antisymmetric");
    // 8 (F^01 F^23 - F^02 F^13 + F^03 F^12) with F^01 = -0.3, F^23 = 0.2.
    assert(abs(scalar(flat, "epsilon_{abcd} F^{ab} F^{cd}") + 8.0 * 0.3 * 0.2) < 1.0e-12, "epsilon_{abcd} F^{ab} F^{cd} is wrong");
    // Contravariant definition stores the same tensor.
    flat := GeneralRelativity.IndexNotation.defineMatrix(flat, "H", eta * F * eta, {true, true});
    assert(maxAbsMatrix(matrix(flat, "H_{ab} - F_{ab}")) < 1.0e-16, "defineMatrix with upper indices did not lower them");
    // Null vector.
    w := {1.0, 0.6, 0.8, 0.0};
    flat := GeneralRelativity.IndexNotation.defineVector(flat, "k", w);
    assert(abs(scalar(flat, "k_a k^a")) < 1.0e-15, "k is not null");
    assert(maxAbsVector(vector(flat, "F^{a}_{b} k^b") - eta * F * w) < 1.0e-14, "F^a_b k^b is wrong on flat space");
    assert(maxAbsVector(vector(flat, "F_{ab} k^b") - F * w) < 1.0e-14, "F_ab k^b is wrong on flat space");
    passed := true;
  end runUserTensors;

  // Run by the compiled model at initialization, with time as an argument: a constant call (a binding, or run()
  // without arguments) is evaluated by the omc front end in its interpreter while instantiating, and with this
  // string-parsing code OpenModelica 1.28 then crashes when it exits (exit code 139 after the tests pass).
  discrete Boolean passed(start = false, fixed = true);
  discrete Boolean passedUser(start = false, fixed = true);
algorithm
  when initial() then
    passed := run(time);
    passedUser := runUserTensors(time);
    assert(passed, "GeneralRelativity index-notation assertions did not complete");
    assert(passedUser, "GeneralRelativity index-notation user-tensor assertions did not complete");
  end when;
  annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 1.0e-8, Interval = 0.1));
end GeneralRelativityIndexTests;
