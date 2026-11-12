within Tests;
model GeneralRelativityTests "Curvature pipeline, stress-energy sources, perturbation theory, geodesics and field solvers"
  function run
    output Boolean passed;
  protected
    constant Real pi = 3.141592653589793;
    constant Real rs = 2.0;
    constant Real step[4] = {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3};
    constant Real eta[4, 4] = [-1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
    GeneralRelativity.Tensors.Curvature flat;
    GeneralRelativity.Tensors.Curvature bh;
    GeneralRelativity.Tensors.Curvature kerr0;
    GeneralRelativity.Tensors.Curvature kerr;
    Real x[4];
    Real GammaExact[4, 4, 4];
    Real GammaOrbit[4, 4, 4];
    Real GammaError;
    Real symmetryError;
    Real T[4, 4];
    Real Td[4, 4];
    Real Tr[4, 4];
    Real gInv[4, 4];
    Boolean accepted;
    Real u[4];
    Real h[4, 4];
    Real hBar[4, 4];
    Real hessian[3, 3];
    Real G[4, 4];
    Real Ric[4, 4];
    Real R[4, 4, 4, 4];
    Real laplacian;
    Real ricciScalar;
    Real hessianFD[3, 3];
    Real gradientFD[3];
    Real xs[3];
    Real xPlus[3];
    Real xMinus[3];
    Real p[3];
    Real pTwo[5];
    Real residual[4, 4];
    Real rNode[51];
    Real rhoNode[51];
    Real phiFE[51];
    Real mass[51];
    Real GM;
    Real a[3];
    Real a4[4];
    Real aLT[3];
    Real rhoGrid[9, 9, 9];
    Real boundary[9, 9, 9];
    Real phiGrid[9, 9, 9];
    Real exactGrid[9, 9, 9];
    Integer sweeps;
    Real sorResidual;
    Boolean converged;
    Real gridError;
    Real kWave;
    Real spacing;
    Real position[3];
  algorithm
    // Flat space: every curvature quantity vanishes.
    flat := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.minkowski, {0.3, 1.0, -2.0, 0.5}, {0.0}, step);
    assert(flat.valid > 0.5, "Minkowski metric was not inverted");
    assert(Tests.Assertions.maxAbsMatrix(flat.ricci) < 1.0e-12 and abs(flat.kretschmann) < 1.0e-12,
      "Minkowski space has nonzero curvature");

    // Schwarzschild: Christoffel symbols against closed form, vacuum, Kretschmann.
    x := {0.0, 5.0 * rs, 1.1, 0.4};
    bh := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.schwarzschild, x, {rs}, step);
    GammaExact := GeneralRelativity.Metrics.schwarzschildChristoffel(x, rs);
    GammaError := 0.0;
    for i in 1:4 loop
      for j in 1:4 loop
        for k in 1:4 loop
          GammaError := max(GammaError, abs(bh.Gamma[i, j, k] - GammaExact[i, j, k]));
        end for;
      end for;
    end for;
    assert(GammaError < 1.0e-5, "Schwarzschild Christoffel symbols disagree with the closed form");
    assert(Tests.Assertions.maxAbsMatrix(bh.ricci) < 1.0e-5, "Schwarzschild Ricci tensor is not zero");
    assert(abs(bh.ricciScalar) < 1.0e-6, "Schwarzschild Ricci scalar is not zero");
    assert(abs(bh.kretschmann / (12.0 * rs^2 / x[2]^6) - 1.0) < 1.0e-4,
      "Schwarzschild Kretschmann scalar is not 12 rs^2 / r^6");
    // Riemann symmetries: antisymmetric in each pair, symmetric under pair exchange.
    symmetryError := 0.0;
    for i in 1:4 loop
      for j in 1:4 loop
        for k in 1:4 loop
          for l in 1:4 loop
            symmetryError := max(symmetryError, abs(bh.riemann[i, j, k, l] + bh.riemann[j, i, k, l]));
            symmetryError := max(symmetryError, abs(bh.riemann[i, j, k, l] + bh.riemann[i, j, l, k]));
            symmetryError := max(symmetryError, abs(bh.riemann[i, j, k, l] - bh.riemann[k, l, i, j]));
          end for;
        end for;
      end for;
    end for;
    assert(symmetryError < 1.0e-9, "Riemann tensor lost its index symmetries");
    // Implied stress-energy of a vacuum solution is zero.
    T := GeneralRelativity.Tensors.stressEnergyFromEinstein(bh.einstein, bh.g, 1.0, 1.0);
    assert(Tests.Assertions.maxAbsMatrix(T) < 1.0e-5, "Vacuum metric implied a nonzero stress-energy tensor");

    // Kerr: a = 0 reproduces Schwarzschild, and a spinning hole is still vacuum.
    kerr0 := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.kerr, x, {rs, 0.0}, step);
    assert(Tests.Assertions.maxAbsMatrix(kerr0.g - bh.g) < 1.0e-12, "Kerr with a = 0 is not Schwarzschild");
    kerr := GeneralRelativity.Metrics.curvatureAt(GeneralRelativity.Metrics.kerr, x, {rs, 0.4 * rs}, step);
    assert(kerr.valid > 0.5 and Tests.Assertions.maxAbsMatrix(kerr.ricci) < 1.0e-5, "Kerr Ricci tensor is not zero");
    assert(kerr.kretschmann > 0.0 and abs(kerr.kretschmann / bh.kretschmann - 1.0) < 0.1,
      "Kerr Kretschmann scalar is not close to Schwarzschild at r = 5 rs");

    // Stress-energy sources in the Minkowski rest frame.
    u := {1.0, 0.0, 0.0, 0.0};
    T := GeneralRelativity.StressEnergy.perfectFluid(eta, u, 3.0, 0.5, 1.0);
    assert(Tests.Assertions.maxAbsMatrix(T - diagonal({3.0, 0.5, 0.5, 0.5})) < 1.0e-12,
      "Perfect fluid at rest is not diag(rho c^2, p, p, p)");
    Td := GeneralRelativity.StressEnergy.dust(eta, u, 3.0);
    assert(Tests.Assertions.maxAbsMatrix(Td - diagonal({3.0, 0.0, 0.0, 0.0})) < 1.0e-12, "Dust at rest is not diag(rho c^2, 0, 0, 0)");
    Tr := GeneralRelativity.StressEnergy.radiation(eta, u, 3.0, 1.0);
    assert(abs(GeneralRelativity.StressEnergy.trace(eta, Tr)) < 1.0e-12, "Radiation stress-energy is not traceless");
    assert(abs(GeneralRelativity.StressEnergy.trace(eta, T) - (-3.0 + 1.5)) < 1.0e-12, "Perfect fluid trace is not -rho c^2 + 3p");
    assert(abs(GeneralRelativity.StressEnergy.energyDensity(T, u, 1.0) - 3.0) < 1.0e-12, "Fluid-frame energy density is not rho");
    // A boosted fluid keeps its rest-frame energy density.
    u := GeneralRelativity.Geodesics.fourVelocity(eta, {0.6, 0.0, 0.0}, 1.0);
    assert(abs(GeneralRelativity.Geodesics.normResidual(eta, u, 1.0)) < 1.0e-12, "fourVelocity is not normalized");
    assert(abs(u[1] - 1.25) < 1.0e-12, "Lorentz factor of v = 0.6 c is not 1.25");
    T := GeneralRelativity.StressEnergy.perfectFluid(eta, u, 3.0, 0.5, 1.0);
    assert(abs(GeneralRelativity.StressEnergy.energyDensity(T, u, 1.0) - 3.0) < 1.0e-12,
      "Boosted perfect fluid changed its rest-frame energy density");
    // A cosmological constant on the left-hand side equals the vacuum-energy source on the right.
    T := GeneralRelativity.StressEnergy.cosmologicalConstant(bh.g, 0.01, 1.0, 1.0);
    residual := GeneralRelativity.StressEnergy.einsteinResidual(bh.einstein, bh.g, T, 1.0, 1.0, 0.0)
      - GeneralRelativity.StressEnergy.einsteinResidual(bh.einstein, bh.g, zeros(4, 4), 1.0, 1.0, 0.01);
    assert(Tests.Assertions.maxAbsMatrix(residual) < 1.0e-12, "Cosmological-constant source is not equivalent to the Lambda term");

    // Perturbation theory: trace reversal is an involution; weak field round trip.
    h := [0.2, 0.1, 0.0, 0.0; 0.1, 0.3, 0.05, 0.0; 0.0, 0.05, -0.1, 0.0; 0.0, 0.0, 0.0, 0.4];
    hBar := GeneralRelativity.Perturbation.traceReverse(GeneralRelativity.Perturbation.traceReverse(h));
    assert(Tests.Assertions.maxAbsMatrix(hBar - h) < 1.0e-12, "traceReverse applied twice is not the identity");
    assert(abs(GeneralRelativity.Perturbation.newtonianPotential(
      GeneralRelativity.Perturbation.weakFieldMetric(-7.0, zeros(3), 3.0), 3.0) + 7.0) < 1.0e-12,
      "newtonianPotential does not invert weakFieldMetric");
    // Linearized Einstein tensor of a static potential: G_00 = 2 laplacian / c^2, G_ij = 0.
    hessian := [1.0, 0.2, -0.3; 0.2, -2.0, 0.4; -0.3, 0.4, 0.5];
    (R, G, laplacian) := GeneralRelativity.Perturbation.weakFieldCurvature(hessian, 2.0);
    assert(abs(G[1, 1] - 2.0 * (-0.5) / 4.0) < 1.0e-12, "Linearized G_00 is not 2 laplacian(Phi) / c^2");
    assert(Tests.Assertions.maxAbsMatrix(G[2:4, 2:4]) < 1.0e-12 and Tests.Assertions.maxAbsVector(G[1, 2:4]) < 1.0e-12,
      "Linearized G_ij and G_0i of a static potential are not zero");
    assert(Tests.Assertions.maxAbsMatrix(GeneralRelativity.Perturbation.tidalTensor(R, 2.0) - hessian) < 1.0e-12,
      "Tidal tensor c^2 R_0i0j is not the Hessian of the potential");
    // Point mass Hessian is traceless and matches the differenced gradient.
    xs := {1.0, -2.0, 0.5};
    hessian := GeneralRelativity.Perturbation.Potentials.pointMassHessian(xs, 3.0);
    assert(abs(hessian[1, 1] + hessian[2, 2] + hessian[3, 3]) < 1.0e-12, "Point-mass tidal tensor is not traceless");
    hessianFD := GeneralRelativity.Perturbation.Potentials.hessianFromGradient(
      GeneralRelativity.Perturbation.Potentials.pointMassGradient, xs, {3.0}, 1.0e-4);
    assert(Tests.Assertions.maxAbsMatrix(hessianFD - hessian) < 1.0e-6, "hessianFromGradient disagrees with pointMassHessian");
    // J2 and two-body gradients are the gradients of their potentials.
    p := {3.0, 0.7, 0.01};
    gradientFD := zeros(3);
    for i in 1:3 loop
      xPlus := xs;
      xMinus := xs;
      xPlus[i] := xs[i] + 1.0e-5;
      xMinus[i] := xs[i] - 1.0e-5;
      gradientFD[i] := (GeneralRelativity.Perturbation.Potentials.oblateJ2(xPlus, p)
        - GeneralRelativity.Perturbation.Potentials.oblateJ2(xMinus, p)) / 2.0e-5;
    end for;
    assert(Tests.Assertions.maxAbsVector(gradientFD - GeneralRelativity.Perturbation.Potentials.oblateJ2Gradient(xs, p)) < 1.0e-7,
      "oblateJ2Gradient is not the gradient of oblateJ2");
    pTwo := {3.0, 0.5, 4.0, 1.0, -1.0};
    for i in 1:3 loop
      xPlus := xs;
      xMinus := xs;
      xPlus[i] := xs[i] + 1.0e-5;
      xMinus[i] := xs[i] - 1.0e-5;
      gradientFD[i] := (GeneralRelativity.Perturbation.Potentials.twoBody(xPlus, pTwo)
        - GeneralRelativity.Perturbation.Potentials.twoBody(xMinus, pTwo)) / 2.0e-5;
    end for;
    assert(Tests.Assertions.maxAbsVector(gradientFD - GeneralRelativity.Perturbation.Potentials.twoBodyGradient(xs, pTwo)) < 1.0e-7,
      "twoBodyGradient is not the gradient of twoBody");
    // The weak-field Earth metric encodes the J2 potential in g_00 and the Lense-Thirring terms in g_0i.
    position := {7.0e6, 1.0e6, 2.0e6};
    h := GeneralRelativity.Metrics.earth({0.0, position[1], position[2], position[3]},
      {GeneralRelativity.Constants.GM_earth, GeneralRelativity.Constants.R_earth, GeneralRelativity.Constants.J2_earth,
       GeneralRelativity.Constants.c, GeneralRelativity.Constants.J_earth});
    assert(abs(GeneralRelativity.Perturbation.newtonianPotential(h, GeneralRelativity.Constants.c)
      - GeneralRelativity.Perturbation.Potentials.oblateJ2(position,
        {GeneralRelativity.Constants.GM_earth, GeneralRelativity.Constants.R_earth, GeneralRelativity.Constants.J2_earth}))
      < 1.0e-5 * GeneralRelativity.Constants.GM_earth / 7.0e6,
      "Metrics.earth g_00 does not encode the J2 potential (1 + g_00 keeps about eight digits of it)");
    // h_0i = k J (y, -x, 0) with k = 2G/(c^3 r^3) > 0: zero z component, orthogonal to (x, y), and h_0x has the sign of y.
    assert(abs(h[1, 4]) < 1.0e-30 and h[1, 2] > 0.0 and h[1, 3] < 0.0
      and abs(h[1, 2] * position[1] + h[1, 3] * position[2]) < 1.0e-9 * abs(h[1, 2] * position[1]),
      "Lense-Thirring g_0i is not -2G(J x r)/(c^3 r^3) for J along z");
    // Post-Newtonian acceleration on a circular orbit is radial: (4GM/r - v^2) x / r^3 * GM / c^2.
    a := GeneralRelativity.Perturbation.postNewtonianAcceleration({1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 1.0, 10.0);
    assert(Tests.Assertions.maxAbsVector(a - {0.03, 0.0, 0.0}) < 1.0e-12, "postNewtonianAcceleration is wrong on a circular orbit");
    aLT := GeneralRelativity.Perturbation.lenseThirringAcceleration({1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, 1.0, 1.0);
    assert(Tests.Assertions.maxAbsVector(aLT - {2.0, 0.0, 0.0}) < 1.0e-12,
      "lenseThirringAcceleration in the equatorial plane is not 2 G J v / (c^2 r^3) radially outward");

    // Geodesics: Killing constants of the Schwarzschild metric and a circular-orbit acceleration.
    (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(bh.g);
    u := {1.5, 0.0, 0.0, 0.02};
    assert(abs(-GeneralRelativity.Geodesics.killingConserved(bh.g, {1.0, 0.0, 0.0, 0.0}, u) - (1.0 - rs / x[2]) * 1.5) < 1.0e-12,
      "Killing energy is not (1 - rs/r) u^t");
    // Circular equatorial orbit: u^phi / u^t = sqrt(M / r^3) makes the radial acceleration vanish.
    x := {0.0, 5.0 * rs, pi / 2.0, 0.0};
    u := {1.0, 0.0, 0.0, sqrt(0.5 * rs / x[2]^3)};
    (a4, h, GammaOrbit) := GeneralRelativity.Geodesics.geodesicAcceleration(GeneralRelativity.Metrics.schwarzschild, x, {rs}, step, u);
    assert(Tests.Assertions.maxAbsVector(a4) < 1.0e-9,
      "Circular Schwarzschild orbit has nonzero geodesic acceleration");

    // Spherical finite elements: density rho = 1 - r for r < 1 (piecewise linear, so interpolated exactly),
    // M = pi/3, Phi(0) = -2 pi/3 = -2 GM, Phi(1) = -GM.
    for i in 1:51 loop
      rNode[i] := 2.0 * (i - 1) / 50.0;
      rhoNode[i] := max(1.0 - rNode[i], 0.0);
    end for;
    (phiFE, mass) := GeneralRelativity.Numerical.poissonFiniteElementSpherical(rNode, rhoNode, 1.0);
    GM := pi / 3.0;
    assert(abs(mass[51] / GM - 1.0) < 1.0e-12, "Finite-element enclosed mass is not pi/3");
    assert(abs(phiFE[1] / (-2.0 * GM) - 1.0) < 5.0e-3, "Finite-element central potential is not -2 GM");
    assert(abs(phiFE[26] / (-GM) - 1.0) < 1.0e-3, "Finite-element surface potential is not -GM/R");
    assert(abs(phiFE[51] / (-0.5 * GM) - 1.0) < 1.0e-12, "Finite-element outer boundary value is wrong");

    // Finite-difference Poisson with a manufactured solution: Phi = sin(kx) sin(ky) sin(kz) on [0, pi/k]^3,
    // rho = -3 k^2 Phi / (4 pi G), zero Dirichlet boundary. The discrete error is O(spacing^2).
    kWave := 1.0;
    spacing := pi / 8.0;
    for i in 1:9 loop
      for j in 1:9 loop
        for l in 1:9 loop
          exactGrid[i, j, l] := sin(kWave * (i - 1) * spacing) * sin(kWave * (j - 1) * spacing) * sin(kWave * (l - 1) * spacing);
          rhoGrid[i, j, l] := -3.0 * kWave^2 * exactGrid[i, j, l] / (4.0 * pi);
          boundary[i, j, l] := 0.0;
        end for;
      end for;
    end for;
    (phiGrid, sweeps, sorResidual, converged) := GeneralRelativity.Numerical.poissonFiniteDifference3D(
      rhoGrid, {spacing, spacing, spacing}, boundary, 1.0, 1.5, 2000, 1.0e-12);
    assert(converged, "SOR did not converge on the manufactured Poisson problem");
    gridError := 0.0;
    for i in 1:9 loop
      for j in 1:9 loop
        for l in 1:9 loop
          gridError := max(gridError, abs(phiGrid[i, j, l] - exactGrid[i, j, l]));
        end for;
      end for;
    end for;
    assert(gridError < 0.03, "Finite-difference Poisson solution misses the manufactured solution");

    passed := true;
  end run;

  parameter Boolean passed = run();
equation
  assert(passed, "GeneralRelativity assertions did not complete");
  annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 1.0e-8, Interval = 0.001));
end GeneralRelativityTests;
