within Tests;
model GeneralRelativityCompactBinaryTests "Compact binaries: post-Newtonian circular orbits and radiation reaction, remnant fits, ringdown, quadrupole strain and the superposed Kerr-Schild metric"
  function run
    output Boolean passed;
  protected
    Real eta;
    Real r;
    Real gamma;
    Real speed;
    Real omega;
    Real rel[3];
    Real a[3];
    Real aRR[3];
    Real h[2];
    Real phi;
    Real gB[4, 4];
    Real gK[4, 4];
    Real gS[4, 4];
    Real err;
    Real X[4];
    Real T1;
    Real T2;
    Real radii[3] = {50.0, 100.0, 200.0};
    Real points[4, 4] = {{0.0, 3.0, 1.0, 2.0}, {1.0, 5.0, -2.0, 0.5}, {2.0, -1.2, 0.7, -0.4}, {0.0, 0.3, 0.2, 9.0}};
  algorithm
    // Circular orbits: the exact circular speed of the 2PN conservative dynamics balances the radial pull,
    // and its frequency agrees with the 2PN Kepler law up to O(gamma^3).
    eta := 0.25;
    for k in 1:3 loop
      r := radii[k];
      gamma := 1.0 / r;
      speed := GeneralRelativity.CompactBinaries.circularOrbitSpeed(r, 1.0, eta, 2);
      a := GeneralRelativity.CompactBinaries.relativeAcceleration({r, 0.0, 0.0}, {0.0, speed, 0.0}, 1.0, eta, 2, false);
      assert(abs(a[1] * r + speed^2) < 1.0e-14 and abs(a[2]) < 1.0e-18, "circularOrbitSpeed does not balance the 2PN radial acceleration");
      omega := speed / r;
      rel[k] := omega^2 / GeneralRelativity.CompactBinaries.circularOmega(r, 1.0, eta, 2)^2 - 1.0;
      assert(abs(rel[k]) < 30.0 * gamma^3, "2PN circular frequency differs from the Kepler law at O(gamma^2) or below");
    end for;
    // The residual is a gamma^3 term: halving gamma divides it by 8.
    assert(abs(rel[2] / rel[3] / 8.0 - 1.0) < 0.05 and abs(rel[1] / rel[2] / 8.0 - 1.0) < 0.1, "2PN circular-frequency residual does not scale like gamma^3");
    // The 1PN law alone: residual O(gamma^2).
    speed := GeneralRelativity.CompactBinaries.circularOrbitSpeed(200.0, 1.0, eta, 1);
    assert(abs((speed / 200.0)^2 / GeneralRelativity.CompactBinaries.circularOmega(200.0, 1.0, eta, 1)^2 - 1.0) < 20.0 / 200.0^2, "1PN circular frequency is wrong");
    // Newtonian Kepler law.
    assert(abs(GeneralRelativity.CompactBinaries.circularOrbitSpeed(7.0, 1.0, eta, 0) - sqrt(1.0 / 7.0)) < 1.0e-15, "Newtonian circular speed is wrong");

    // Radiation reaction on a Newtonian circular orbit removes energy at the quadrupole rate:
    // v . a_RR = -(32/5) eta M^4 / r^5 per unit reduced mass, i.e. dE/dt = -(32/5) eta^2 M^5 / r^5.
    for k in 1:2 loop
      eta := if k == 1 then 0.25 else 0.16;
      r := 20.0;
      speed := sqrt(1.0 / r);
      phi := 0.7;
      a := GeneralRelativity.CompactBinaries.relativeAcceleration({r * cos(phi), r * sin(phi), 0.0}, {-speed * sin(phi), speed * cos(phi), 0.0}, 1.0, eta, 0, true);
      aRR := a - GeneralRelativity.CompactBinaries.relativeAcceleration({r * cos(phi), r * sin(phi), 0.0}, {-speed * sin(phi), speed * cos(phi), 0.0}, 1.0, eta, 0, false);
      assert(abs((aRR * {-speed * sin(phi), speed * cos(phi), 0.0}) / (-32.0 / 5.0 * eta / r^5) - 1.0) < 1.0e-12, "2.5PN radiation reaction does not give the quadrupole energy flux");
    end for;
    // The leading-order inspiral: Peters' coalescence time integrates the radial drift, dT/dr0 = -1/rdot(r0).
    T1 := GeneralRelativity.CompactBinaries.mergerTimePeters(10.0 + 1.0e-4, 1.0, 0.2);
    T2 := GeneralRelativity.CompactBinaries.mergerTimePeters(10.0 - 1.0e-4, 1.0, 0.2);
    assert(abs((T1 - T2) / 2.0e-4 * GeneralRelativity.CompactBinaries.inspiralRadialVelocity(10.0, 1.0, 0.2) + 1.0) < 1.0e-8, "Peters time and radial drift are inconsistent");

    // Remnant of an equal-mass merger of non-spinning holes; the test-mass limit.
    assert(abs(GeneralRelativity.CompactBinaries.finalSpin(0.25) - 0.6865) < 5.0e-4, "Equal-mass final spin is not 0.686");
    assert(abs(GeneralRelativity.CompactBinaries.radiatedEnergyFraction(0.25) - 0.0483) < 5.0e-4, "Equal-mass radiated energy is not 4.8%");
    assert(abs(GeneralRelativity.CompactBinaries.radiatedEnergyFraction(1.0e-6) / 1.0e-6 - (1.0 - sqrt(8.0 / 9.0))) < 1.0e-5,
      "Small mass ratio: the radiated energy is not the binding energy at the innermost stable orbit");
    // Ringdown of a Schwarzschild hole: M omega = 0.37367, tau = 11.24 M (Leaver); the fits are within a few percent.
    assert(abs(GeneralRelativity.CompactBinaries.ringdownOmega(0.0) / 0.37367 - 1.0) < 0.02, "Schwarzschild ringdown frequency is off");
    assert(abs(2.0 * GeneralRelativity.CompactBinaries.ringdownQuality(0.0) / GeneralRelativity.CompactBinaries.ringdownOmega(0.0) / 11.24 - 1.0) < 0.03,
      "Schwarzschild ringdown damping time is off");
    assert(GeneralRelativity.CompactBinaries.ringdownOmega(0.9) > GeneralRelativity.CompactBinaries.ringdownOmega(0.5)
      and GeneralRelativity.CompactBinaries.ringdownQuality(0.9) > GeneralRelativity.CompactBinaries.ringdownQuality(0.5), "Faster spin should ring higher and longer");

    // Quadrupole strain of a Newtonian circular orbit: amplitude 4 eta (M omega)^(2/3) at any phase.
    eta := 0.2;
    r := 10.0;
    omega := sqrt(1.0 / r^3);
    for k in 0:3 loop
      phi := 0.4 + 1.1 * k;
      h := GeneralRelativity.CompactBinaries.quadrupoleStrain({r * cos(phi), r * sin(phi), 0.0}, {-r * omega * sin(phi), r * omega * cos(phi), 0.0}, 1.0, eta);
      assert(abs(sqrt(h * h) / (4.0 * eta * omega^(2.0 / 3.0)) - 1.0) < 1.0e-9, "Quadrupole strain amplitude is wrong");
      // h+ - i hx = -A exp(-2 i phi): the wave turns at twice the orbital rate.
      assert(abs(h[1] + 4.0 * eta * omega^(2.0 / 3.0) * cos(2.0 * phi)) < 1.0e-12 and abs(h[2] + 4.0 * eta * omega^(2.0 / 3.0) * sin(2.0 * phi)) < 1.0e-12,
        "Quadrupole strain phase is wrong");
    end for;

    // Superposed Kerr-Schild: one body alone is Schwarzschild in Kerr-Schild form (rs = 2 m), anywhere;
    // the metric is symmetric under exchanging the bodies.
    for k in 1:4 loop
      X := points[k, :];
      gB := GeneralRelativity.Metrics.binaryKerrSchild(X, {0.7, 0.0, 0.0, 0.0, 0.0, 20.0, -3.0, 1.0});
      gK := GeneralRelativity.Metrics.kerrSchild(X, {1.4});
      err := max(abs(gB - gK));
      assert(err < 1.0e-12, "binaryKerrSchild with one body is not kerrSchild");
      gB := GeneralRelativity.Metrics.binaryKerrSchild(X + {0.0, 1.0, -2.0, 0.5}, {0.0, 0.7, 6.0, 6.0, 6.0, 1.0, -2.0, 0.5});
      assert(max(abs(gB - gK)) < 1.0e-12, "binaryKerrSchild is not translation-covariant");
      gB := GeneralRelativity.Metrics.binaryKerrSchild(X, {0.6, 0.4, 1.5, 0.2, 0.0, -2.0, -0.3, 0.1});
      gS := GeneralRelativity.Metrics.binaryKerrSchild(X, {0.4, 0.6, -2.0, -0.3, 0.1, 1.5, 0.2, 0.0});
      assert(max(abs(gB - gS)) < 1.0e-15, "binaryKerrSchild is not symmetric under exchanging the bodies");
    end for;
    passed := true;
  end run;

  parameter Boolean passed = run();
equation
  assert(passed, "Compact-binary assertions did not complete");
  annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 1.0e-8, Interval = 0.001));
end GeneralRelativityCompactBinaryTests;
