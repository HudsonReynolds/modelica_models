within Tests;
model GeneralRelativityConnectionTests "Dimension-generic curvature, general affine connections (torsion, contortion, non-metricity) and Newton-Cartan gravity"
  function schwarzschildConnection "Test connection: closed-form Schwarzschild Christoffel symbols, p = {rs}"
    extends GeneralRelativity.Interfaces.ConnectionFunction;
  algorithm
    Gamma := GeneralRelativity.Metrics.schwarzschildChristoffel(x, p[1]);
  end schwarzschildConnection;

  function uniformSphereConnection "Test connection: Newton-Cartan connection inside a uniform ball, Phi = (2 pi / 3) G rho r^2; p = {rho}, G = c = 1"
    extends GeneralRelativity.Interfaces.ConnectionFunction;
  algorithm
    Gamma := zeros(4, 4, 4);
    for i in 1:3 loop
      Gamma[i + 1, 1, 1] := 4.0 * GeneralRelativity.Constants.pi / 3.0 * p[1] * x[i + 1];
    end for;
  end uniformSphereConnection;

  function weylConnection "Test connection: Minkowski Levi-Civita plus the Weyl term 1/2 (delta^i_j A_k + delta^i_k A_j - eta_jk A^i); p = A_a"
    extends GeneralRelativity.Interfaces.ConnectionFunction;
  protected
    Real eta[4, 4] = diagonal({-1.0, 1.0, 1.0, 1.0});
    Real AUp[4] = eta * p[1:4];
  algorithm
    for i in 1:4 loop
      for j in 1:4 loop
        for k in 1:4 loop
          Gamma[i, j, k] := 0.5 * ((if i == j then p[k] else 0.0) + (if i == k then p[j] else 0.0) - eta[j, k] * AUp[i]);
        end for;
      end for;
    end for;
  end weylConnection;

  function run
    output Boolean passed;
  protected
    constant Real pi = GeneralRelativity.Constants.pi;
    constant Real rs = 2.0;
    constant Real step[4] = {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3};
    Real a;
    Real theta;
    Real g2[2, 2];
    Real gInv2[2, 2];
    Real dg2[2, 2, 2];
    Real ddg2[2, 2, 2, 2];
    Real Gamma2[2, 2, 2];
    Real R2[2, 2, 2, 2];
    Real Ric2[2, 2];
    Real scalar2;
    Real g3[3, 3];
    Real gInv3[3, 3];
    Real dg3[3, 3, 3];
    Real ddg3[3, 3, 3, 3];
    Real Gamma3[3, 3, 3];
    Real R3[3, 3, 3, 3];
    Real r;
    Boolean accepted;
    Real x[4];
    GeneralRelativity.Tensors.Curvature bh;
    GeneralRelativity.Connections.ConnectionCurvature kc;
    Real Gamma[4, 4, 4];
    Real GammaExact[4, 4, 4];
    Real RMixed[4, 4, 4, 4];
    Real dg[4, 4, 4];
    Real Q[4, 4, 4];
    Real T[4, 4, 4];
    Real T2[4, 4, 4];
    Real K[4, 4, 4];
    Real KLower[4, 4, 4];
    Real eta[4, 4];
    Real A[4];
    Real clockResidual;
    Real cometricResidual;
    Real h[4, 4];
    Real tidal[3, 3];
    Real position[3];
    Real residual[4, 4];
    Real u[4];
    Real acc[4];
    Real GammaNC[4, 4, 4];
    Real relative;
  algorithm
    // Dimension-generic Tensors: a 2-sphere of radius a with analytic derivatives, coordinates {theta, phi}.
    a := 3.0;
    theta := 0.7;
    g2 := diagonal({a^2, a^2 * sin(theta)^2});
    dg2 := zeros(2, 2, 2);
    dg2[1, 2, 2] := 2.0 * a^2 * sin(theta) * cos(theta);
    ddg2 := zeros(2, 2, 2, 2);
    ddg2[1, 1, 2, 2] := 2.0 * a^2 * cos(2.0 * theta);
    (gInv2, accepted) := GeneralRelativity.Tensors.inverseMetric(g2);
    Gamma2 := GeneralRelativity.Tensors.christoffel(gInv2, dg2);
    assert(abs(Gamma2[1, 2, 2] + sin(theta) * cos(theta)) < 1.0e-12 and abs(Gamma2[2, 1, 2] - cos(theta) / sin(theta)) < 1.0e-12,
      "2-sphere Christoffel symbols are wrong");
    R2 := GeneralRelativity.Tensors.riemannLower(g2, ddg2, Gamma2);
    assert(abs(R2[1, 2, 1, 2] - a^2 * sin(theta)^2) < 1.0e-12, "2-sphere R_thphthph is not a^2 sin^2 theta");
    Ric2 := GeneralRelativity.Tensors.ricci(gInv2, R2);
    scalar2 := GeneralRelativity.Tensors.ricciScalar(gInv2, Ric2);
    assert(abs(scalar2 - 2.0 / a^2) < 1.0e-12, "2-sphere Ricci scalar is not 2 / a^2");
    assert(abs(GeneralRelativity.Tensors.kretschmann(gInv2, R2) - 4.0 / a^4) < 1.0e-12, "2-sphere Kretschmann scalar is not R^2 = 4 / a^4");
    assert(Tests.Assertions.maxAbsMatrix(GeneralRelativity.Tensors.einstein(g2, Ric2, scalar2)) < 1.0e-12,
      "The Einstein tensor of any 2-dimensional metric vanishes");
    // Flat 3-space in spherical coordinates {r, theta, phi}: all curvature is zero.
    r := 2.5;
    g3 := diagonal({1.0, r^2, r^2 * sin(theta)^2});
    dg3 := zeros(3, 3, 3);
    dg3[1, 2, 2] := 2.0 * r;
    dg3[1, 3, 3] := 2.0 * r * sin(theta)^2;
    dg3[2, 3, 3] := 2.0 * r^2 * sin(theta) * cos(theta);
    ddg3 := zeros(3, 3, 3, 3);
    ddg3[1, 1, 2, 2] := 2.0;
    ddg3[1, 1, 3, 3] := 2.0 * sin(theta)^2;
    ddg3[1, 2, 3, 3] := 4.0 * r * sin(theta) * cos(theta);
    ddg3[2, 1, 3, 3] := ddg3[1, 2, 3, 3];
    ddg3[2, 2, 3, 3] := 2.0 * r^2 * cos(2.0 * theta);
    (gInv3, accepted) := GeneralRelativity.Tensors.inverseMetric(g3);
    Gamma3 := GeneralRelativity.Tensors.christoffel(gInv3, dg3);
    R3 := GeneralRelativity.Tensors.riemannLower(g3, ddg3, Gamma3);
    assert(abs(GeneralRelativity.Tensors.kretschmann(gInv3, R3)) < 1.0e-20 and
      Tests.Assertions.maxAbsMatrix(GeneralRelativity.Tensors.ricci(gInv3, R3)) < 1.0e-12,
      "Flat 3-space in spherical coordinates has curvature");

    // Connections: the Levi-Civita connection of Schwarzschild, its curvature from Connections.riemann,
    // agrees with the metric pipeline; no torsion, no non-metricity.
    x := {0.0, 5.0 * rs, 1.1, 0.4};
    Gamma := GeneralRelativity.Connections.leviCivita(GeneralRelativity.Metrics.schwarzschild, x, {rs}, step);
    GammaExact := GeneralRelativity.Metrics.schwarzschildChristoffel(x, rs);
    assert(max(abs(Gamma[i, j, k] - GammaExact[i, j, k]) for i in 1:4, j in 1:4, k in 1:4) < 1.0e-5,
      "Connections.leviCivita disagrees with the closed-form Schwarzschild Christoffel symbols");
    bh := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.schwarzschild, x, {rs}, step);
    RMixed := GeneralRelativity.Tensors.riemannMixed(bh.gInv, bh.riemann);
    kc := GeneralRelativity.Connections.curvatureAt(schwarzschildConnection, x, {rs}, step);
    // Both sides are finite differences (step 1e-3: nested first differences against second differences of the
    // metric); they agree to about 2e-5 of the largest component.
    assert(max(abs(kc.riemann[i, j, k, l] - RMixed[i, j, k, l]) for i in 1:4, j in 1:4, k in 1:4, l in 1:4)
      < 5.0e-5 * max(abs(RMixed[i, j, k, l]) for i in 1:4, j in 1:4, k in 1:4, l in 1:4),
      "Curvature of the Schwarzschild connection disagrees with the metric Riemann tensor");
    assert(Tests.Assertions.maxAbsMatrix(kc.ricci) < 1.0e-6, "Ricci tensor of the Schwarzschild connection is not zero");
    assert(max(abs(kc.torsion[i, j, k]) for i in 1:4, j in 1:4, k in 1:4) < 1.0e-15, "Levi-Civita connection has torsion");
    dg := GeneralRelativity.Metrics.derivative(GeneralRelativity.Metrics.schwarzschild, x, {rs}, step);
    Q := GeneralRelativity.Connections.nonMetricity(bh.g, dg, GammaExact);
    // dg is a central difference (relative error about h^2 = 1e-6); Gamma is exact.
    assert(max(abs(Q[i, j, k]) for i in 1:4, j in 1:4, k in 1:4) < 2.0e-6 * max(abs(dg[i, j, k]) for i in 1:4, j in 1:4, k in 1:4),
      "Levi-Civita connection is not metric-compatible");

    // Contortion: adding K(T) to the Christoffel symbols gives a metric connection whose torsion is T.
    T := zeros(4, 4, 4);
    for i in 1:4 loop
      for j in 1:4 loop
        for k in j + 1:4 loop
          T[i, j, k] := 0.1 * sin(1.0 + i + 2.0 * j + 3.0 * k);
          T[i, k, j] := -T[i, j, k];
        end for;
      end for;
    end for;
    K := GeneralRelativity.Connections.contortion(bh.g, bh.gInv, T);
    Gamma := GammaExact + K;
    T2 := GeneralRelativity.Connections.torsion(Gamma);
    assert(max(abs(T2[i, j, k] - T[i, j, k]) for i in 1:4, j in 1:4, k in 1:4) < 1.0e-12,
      "Christoffel plus contortion does not have the prescribed torsion");
    Q := GeneralRelativity.Connections.nonMetricity(bh.g, dg, Gamma);
    assert(max(abs(Q[i, j, k]) for i in 1:4, j in 1:4, k in 1:4) < 2.0e-6 * max(abs(dg[i, j, k]) for i in 1:4, j in 1:4, k in 1:4),
      "Christoffel plus contortion is not metric-compatible");
    for i in 1:4 loop
      for j in 1:4 loop
        for k in 1:4 loop
          KLower[i, j, k] := sum(bh.g[i, l] * K[l, j, k] for l in 1:4);
        end for;
      end for;
    end for;
    assert(max(abs(KLower[i, j, k] + KLower[k, j, i]) for i in 1:4, j in 1:4, k in 1:4) < 1.0e-12,
      "Lowered contortion is not antisymmetric in its outer indices");

    // Weyl geometry: Q_lmn = -A_l g_mn for the Weyl connection, and it is torsion-free.
    eta := diagonal({-1.0, 1.0, 1.0, 1.0});
    A := {0.3, -0.2, 0.5, 0.1};
    Gamma := weylConnection(x, A);
    Q := GeneralRelativity.Connections.nonMetricity(eta, zeros(4, 4, 4), Gamma);
    assert(max(abs(Q[l, m, n] + A[l] * eta[m, n]) for l in 1:4, m in 1:4, n in 1:4) < 1.0e-12,
      "Weyl connection non-metricity is not -A_l g_mn");
    T2 := GeneralRelativity.Connections.torsion(Gamma);
    assert(max(abs(T2[i, j, k]) for i in 1:4, j in 1:4, k in 1:4) < 1.0e-15,
      "Weyl connection has torsion");

    // Newton-Cartan: compatibility with the clock form and the spatial cometric.
    h := GeneralRelativity.NewtonCartan.spatialCometric(x, true);
    assert(Tests.Assertions.maxAbsVector(h * GeneralRelativity.NewtonCartan.clockForm()) < 1.0e-15,
      "Spatial cometric is not degenerate along the clock form");
    assert(Tests.Assertions.maxAbsMatrix(h[2:4, 2:4] * GeneralRelativity.NewtonCartan.spatialMetric(x, true) - identity(3)) < 1.0e-12,
      "spatialMetric is not the inverse of the spatial cometric");
    (clockResidual, cometricResidual) := GeneralRelativity.NewtonCartan.compatibilityAt(
      GeneralRelativity.NewtonCartan.pointMassSpherical, x, {rs}, step, true);
    assert(clockResidual < 1.0e-15 and cometricResidual < 1.0e-6, "Spherical Newton-Cartan connection is not compatible with tau and h");
    x := {0.0, 3.0, -4.0, 2.0};
    (clockResidual, cometricResidual) := GeneralRelativity.NewtonCartan.compatibilityAt(
      GeneralRelativity.NewtonCartan.pointMassCartesian, x, {rs, 0.01, -0.02, 0.03}, step, false);
    assert(clockResidual < 1.0e-15 and cometricResidual < 1.0e-15, "Rotating Newton-Cartan connection is not compatible with tau and h");
    // A rotating chart without mass is flat: Coriolis and centrifugal terms are not curvature.
    kc := GeneralRelativity.Connections.curvatureAt(GeneralRelativity.NewtonCartan.pointMassCartesian, x, {0.0, 0.01, -0.02, 0.03}, step);
    assert(max(abs(kc.riemann[i, j, k, l]) for i in 1:4, j in 1:4, k in 1:4, l in 1:4) < 1.0e-12,
      "Newton-Cartan connection of an empty rotating chart is curved");
    assert(max(abs(kc.torsion[i, j, k]) for i in 1:4, j in 1:4, k in 1:4) < 1.0e-15, "Newton-Cartan connection has torsion");
    // With mass, R^i_0j0 is the tidal tensor Phi_,ij / c^2 whatever the rotation; Ricci vanishes in vacuum.
    kc := GeneralRelativity.Connections.curvatureAt(GeneralRelativity.NewtonCartan.pointMassCartesian, x, {rs, 0.01, -0.02, 0.03}, step);
    position := x[2:4];
    tidal := GeneralRelativity.Perturbation.Potentials.pointMassHessian(position, rs / 2.0);
    assert(max(abs(kc.riemann[i + 1, 1, j + 1, 1] - tidal[i, j]) for i in 1:3, j in 1:3) < 1.0e-8,
      "Newton-Cartan curvature R^i_0j0 is not the tidal tensor");
    residual := GeneralRelativity.NewtonCartan.fieldEquationResidual(kc.ricci, 0.0);
    assert(Tests.Assertions.maxAbsMatrix(residual) < 1.0e-8, "Newton-Cartan Ricci tensor is not zero in vacuum");
    // Inside matter the field equation R_ab = 4 pi G rho tau_a tau_b is Poisson's equation.
    kc := GeneralRelativity.Connections.curvatureAt(uniformSphereConnection, {0.0, 0.2, 0.1, -0.3}, {0.05}, step);
    assert(Tests.Assertions.maxAbsMatrix(GeneralRelativity.NewtonCartan.fieldEquationResidual(kc.ricci, 0.05)) < 1.0e-10,
      "Uniform-density Newton-Cartan connection violates R_00 = 4 pi G rho");
    residual := GeneralRelativity.NewtonCartan.fieldEquationResidual(kc.ricci, 0.0);
    assert(abs(residual[1, 1] - 4.0 * pi * 0.05) < 1.0e-10,
      "fieldEquationResidual does not detect a missing source");
    // Free fall from rest: d^2 r / d(ct)^2 = -GM / (c^2 r^2) = -rs / (2 r^2).
    x := {0.0, 7.0, pi / 2.0, 0.0};
    GammaNC := GeneralRelativity.NewtonCartan.pointMassSpherical(x, {rs});
    u := {1.0, 0.0, 0.0, 0.0};
    acc := GeneralRelativity.Geodesics.acceleration(GammaNC, u);
    assert(Tests.Assertions.maxAbsVector(acc - {0.0, -rs / (2.0 * 7.0^2), 0.0, 0.0}) < 1.0e-15,
      "Newton-Cartan free fall is not -GM/r^2");
    // Newtonian limit: Schwarzschild Gamma^r_00 = (rs / 2 r^2)(1 - rs/r), so the relative difference is exactly rs/r.
    x := {0.0, 1.0e4 * rs, 1.0, 0.0};
    GammaExact := GeneralRelativity.Metrics.schwarzschildChristoffel(x, rs);
    GammaNC := GeneralRelativity.NewtonCartan.pointMassSpherical(x, {rs});
    relative := (GammaNC[2, 1, 1] - GammaExact[2, 1, 1]) / GammaNC[2, 1, 1];
    assert(abs(relative - rs / x[2]) < 1.0e-12, "Schwarzschild and Newton-Cartan Gamma^r_00 do not differ by rs/r");
    assert(max(abs(GammaExact[i, j, k] - GammaNC[i, j, k]) / max(abs(GammaNC[i, j, k]), 1.0 / x[2]) for i in 1:4, j in 1:4, k in 1:4)
      < 2.0 * rs / x[2], "Schwarzschild connection does not approach Newton-Cartan far from the mass");
    passed := true;
  end run;

  parameter Boolean passed = run();
equation
  assert(passed, "GeneralRelativity connection assertions did not complete");
  annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 1.0e-8, Interval = 0.001),
    Documentation(info="<html>
<p>Checks the dimension-generic curvature chain on a 2-sphere and flat 3-space,
the general-connection functions of
<a href=\"modelica://GeneralRelativity.Connections\">Connections</a> against the
Levi-Civita case and against torsion and non-metricity built by hand, and the
Newton-Cartan connection: compatibility with the clock form and spatial
cometric, flatness of a rotating empty chart, tidal curvature, Poisson's
equation as the field equation, and the Newtonian limit of Schwarzschild.</p>
</html>"));
end GeneralRelativityConnectionTests;
