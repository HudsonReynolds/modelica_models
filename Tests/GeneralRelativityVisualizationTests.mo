within Tests;
model GeneralRelativityVisualizationTests "Physics behind the viewer styles: lapse, river speed, frame dragging and light cones, and the Kerr-Schild metric the ray tracer integrates"
  function run
    output Boolean passed;
  protected
    constant Real pi = GeneralRelativity.Constants.pi;
    constant Real rs = 2.0;
    constant Real step[4] = {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3};
    GeneralRelativity.Visualization.SpacetimePoint s;
    GeneralRelativity.Visualization.SpacetimePoint s2;
    GeneralRelativity.Tensors.Curvature kBL;
    GeneralRelativity.Tensors.Curvature kKS;
    Real r;
    Real f;
    Real a;
    Real rho;
    Real arealRadius;
    Real omega;
    Real Delta;
    Real A;
    Real theta;
    Real g[4, 4];
    Real gInv[4, 4];
    Real eta[4, 4];
    Real k[4];
    Boolean accepted;
  algorithm
    // Schwarzschild, spherical chart: alpha = sqrt(1 - rs/r), river = sqrt(rs/r) (Newtonian escape speed),
    // radial light speed 1 - rs/r, tangential sqrt(1 - rs/r), no frame dragging.
    r := 7.0;
    f := 1.0 - rs / r;
    s := GeneralRelativity.Visualization.spacetimePointAt(GeneralRelativity.Metrics.schwarzschild, {0.0, r, pi / 2.0, 0.3}, {rs}, true);
    assert(s.valid > 0.5 and s.coneValid > 0.5, "Schwarzschild exterior point marked invalid");
    assert(abs(s.lapse - sqrt(f)) < 1.0e-14 and abs(s.timeLag - (1.0 - sqrt(f))) < 1.0e-14, "Schwarzschild lapse is not sqrt(1 - rs/r)");
    assert(abs(s.riverSpeed - sqrt(rs / r)) < 1.0e-14, "Schwarzschild river speed is not sqrt(rs/r)");
    assert(abs(s.dragSpeed) < 1.0e-15, "Schwarzschild has frame dragging");
    assert(abs(s.coneOutward - f) < 1.0e-14 and abs(s.coneInward + f) < 1.0e-14, "Schwarzschild radial light speed is not +-(1 - rs/r)");
    assert(abs(s.conePrograde - sqrt(f)) < 1.0e-14 and abs(s.coneRetrograde + sqrt(f)) < 1.0e-14,
      "Schwarzschild tangential light speed is not +-sqrt(1 - rs/r)");
    // The river reaches c at the horizon and exceeds it inside, where no static observers exist.
    // (Just outside r = rs, where g_rr of the Schwarzschild chart is finite.)
    s := GeneralRelativity.Visualization.spacetimePointAt(GeneralRelativity.Metrics.schwarzschild, {0.0, rs * (1.0 + 1.0e-10), pi / 2.0, 0.0}, {rs}, true);
    assert(abs(s.riverSpeed - 1.0) < 1.0e-9, "River speed is not c at the horizon");
    s := GeneralRelativity.Visualization.spacetimePointAt(GeneralRelativity.Metrics.schwarzschild, {0.0, 0.75 * rs, pi / 2.0, 0.0}, {rs}, true);
    assert(s.valid < 0.5 and s.coneValid < 0.5 and abs(s.riverSpeed - sqrt(4.0 / 3.0)) < 1.0e-14,
      "Inside the horizon the lapse must be invalid and the river faster than light");

    // The same spacetime in isotropic Cartesian coordinates: lapse and river are properties of the static
    // observers, not of the chart, so they match at the same areal radius; isotropy makes the cones round.
    rho := 4.0;
    arealRadius := rho * (1.0 + rs / (4.0 * rho))^2;
    s := GeneralRelativity.Visualization.spacetimePointAt(GeneralRelativity.Metrics.isotropicSchwarzschild,
      {0.0, rho * cos(0.7), rho * sin(0.7), 0.0}, {rs}, false);
    assert(abs(s.lapse - sqrt(1.0 - rs / arealRadius)) < 1.0e-14, "Isotropic-chart lapse differs from the Schwarzschild lapse");
    assert(abs(s.riverSpeed - sqrt(rs / arealRadius)) < 1.0e-14, "Isotropic-chart river speed differs from sqrt(rs/r)");
    assert(abs(s.coneOutward - s.conePrograde) < 1.0e-14 and abs(s.coneInward - s.coneRetrograde) < 1.0e-14,
      "Light cones of a conformally flat slice are not round");

    // Minkowski in Cartesian coordinates: unit lapse, no river, cones at +-c.
    s := GeneralRelativity.Visualization.spacetimePointAt(GeneralRelativity.Metrics.minkowski, {0.0, 1.0, -2.0, 0.0}, {0.0}, false);
    assert(abs(s.lapse - 1.0) < 1.0e-15 and abs(s.timeLag) < 1.0e-15 and abs(s.riverSpeed) < 1.0e-15,
      "Flat space has a time lag or a river");
    assert(abs(s.coneOutward - 1.0) < 1.0e-15 and abs(s.coneRetrograde + 1.0) < 1.0e-15, "Flat-space light speed is not c");

    // Weak field: the time lag keeps its digits although 1 - alpha is about 1e-9.
    s := GeneralRelativity.Visualization.spacetimePointAt(GeneralRelativity.Metrics.schwarzschild, {0.0, 1.0e9 * rs, pi / 2.0, 0.0}, {rs}, true);
    assert(abs(s.timeLag / (0.5e-9) - 1.0) < 1.0e-6, "Weak-field time lag lost precision (expected rs/(2r))");

    // Kerr, equatorial plane in Boyer-Lindquist coordinates: closed-form lapse and frame dragging.
    a := 0.6;
    r := 6.0;
    Delta := r^2 - rs * r + a^2;
    A := (r^2 + a^2)^2 - a^2 * Delta;
    omega := rs * a * r / A;
    s := GeneralRelativity.Visualization.spacetimePointAt(GeneralRelativity.Metrics.kerr, {0.0, r, pi / 2.0, 0.0}, {rs, a}, true);
    assert(abs(s.lapse - sqrt(r^2 * Delta / A)) < 1.0e-14, "Kerr equatorial lapse is not sqrt(r^2 Delta / A)");
    assert(abs(s.dragSpeed - omega * r) < 1.0e-14, "Kerr frame dragging is not omega = rs a r / A");
    assert(s.conePrograde > -s.coneRetrograde, "Frame dragging does not tilt the Kerr light cone prograde");
    // The river reaches c at the static limit (r = rs on the equator), outside the horizon.
    s := GeneralRelativity.Visualization.spacetimePointAt(GeneralRelativity.Metrics.kerr, {0.0, rs, pi / 2.0, 0.0}, {rs, a}, true);
    assert(abs(s.riverSpeed - 1.0) < 1.0e-14 and s.valid > 0.5, "Kerr river is not c at the static limit");

    // Kerr-Schild form (the ray tracer's metric): inverse is eta - f k k, frame dragging has the same sense and size
    // as Boyer-Lindquist, curvature invariants agree.
    eta := diagonal({-1.0, 1.0, 1.0, 1.0});
    g := GeneralRelativity.Metrics.kerrSchild({0.0, 3.0, -1.0, 2.0}, {rs, a, 0.4});
    (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(g);
    assert(Tests.Assertions.maxAbsMatrix(g * gInv - identity(4)) < 1.0e-12, "Kerr-Schild metric inversion failed");
    // At (sqrt(r^2 + a^2), 0, 0) in the equatorial plane the Boyer-Lindquist radius is r.
    s2 := GeneralRelativity.Visualization.spacetimePointAt(GeneralRelativity.Metrics.kerrSchild, {0.0, sqrt(r^2 + a^2), 0.0, 0.0}, {rs, a}, false);
    assert(abs(s2.dragSpeed / sqrt(r^2 + a^2) - omega) < 1.0e-12, "Kerr-Schild frame dragging differs from Boyer-Lindquist");
    assert(abs(s2.riverSpeed - sqrt(rs / r)) < 1.0e-12, "Kerr-Schild g_00 differs from Boyer-Lindquist on the equator");
    // A generic point: (r, theta) in Boyer-Lindquist sits at x = r sin(theta), y = a sin(theta), z = r cos(theta).
    theta := 1.1;
    kBL := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.kerr, {0.0, r, theta, 0.0}, {rs, a}, step);
    kKS := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.kerrSchild,
      {0.0, r * sin(theta), a * sin(theta), r * cos(theta)}, {rs, a}, step);
    assert(Tests.Assertions.maxAbsMatrix(kKS.ricci) < 1.0e-6, "Kerr-Schild Kerr is not vacuum");
    assert(abs(kKS.kretschmann / kBL.kretschmann - 1.0) < 1.0e-4, "Kerr-Schild and Boyer-Lindquist Kretschmann scalars differ");
    // Charged, non-rotating: Reissner-Nordstrom, traceless electromagnetic source, same Kretschmann scalar.
    kBL := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.reissnerNordstrom, {0.0, r, theta, 0.0}, {rs, 0.7}, step);
    kKS := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.kerrSchild,
      {0.0, r * sin(theta), 0.0, r * cos(theta)}, {rs, 0.0, 0.7}, step);
    assert(abs(kKS.ricciScalar) < 1.0e-6, "Kerr-Schild Reissner-Nordstrom has a nonzero Ricci scalar");
    assert(abs(kKS.kretschmann / kBL.kretschmann - 1.0) < 1.0e-4, "Kerr-Schild and Reissner-Nordstrom Kretschmann scalars differ");
    // The Kerr-Schild covector is null for eta, which is what makes the inverse exact.
    r := sqrt(0.5 * (14.0 - a^2 + sqrt((14.0 - a^2)^2 + 4.0 * a^2 * 4.0)));
    k := {1.0, (r * 3.0 + a * (-1.0)) / (r^2 + a^2), (r * (-1.0) - a * 3.0) / (r^2 + a^2), 2.0 / r};
    assert(abs(k * (eta * k)) < 1.0e-14, "Kerr-Schild covector is not null");

    // Embedding: for Schwarzschild the circumferential radius is r (Flamm), checked through g_phiphi.
    kBL := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.schwarzschild, {0.0, 7.0, pi / 2.0, 0.0}, {rs}, step);
    assert(abs(sqrt(kBL.g[4, 4]) - 7.0) < 1.0e-14, "Equatorial circumferential radius of Schwarzschild is not r");
    passed := true;
  end run;

  parameter Boolean passed = run();
equation
  assert(passed, "GeneralRelativity visualization assertions did not complete");
  annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 1.0e-8, Interval = 0.001),
    Documentation(info="<html>
<p>Checks <a href=\"modelica://GeneralRelativity.Visualization.spacetimePointAt\">spacetimePointAt</a>
against closed forms for Schwarzschild (spherical and isotropic charts), Minkowski, a weak field and
Kerr, and the <a href=\"modelica://GeneralRelativity.Metrics.kerrSchild\">Kerr-Schild metric</a> used by the
viewer's ray tracer against the Boyer-Lindquist Kerr and Reissner-Nordstr&ouml;m metrics.</p>
</html>"));
end GeneralRelativityVisualizationTests;
