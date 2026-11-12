within GeneralRelativity.Examples;
model WeakFieldPoisson "Finite-difference and finite-element potentials of a homogeneous sphere, compared with the closed form and lifted to the weak-field metric"
  parameter Real rho = 1.0 "Density in geometric units, G = 1";
  parameter Real radius = 1.0 "Sphere radius";
  parameter Integer n = 61 "Finite-difference nodes per axis";
  parameter Real halfWidth = 4.0 "Half-width of the finite-difference box";
  parameter Integer nElements = 100 "Finite elements on [0, 4 radius]";
  parameter Real c = 1.0e3 "Speed of light in these units, so Phi/c^2 is small";
  final parameter Real GM = 4.0 / 3.0 * GeneralRelativity.Constants.pi * rho * radius^3 "Enclosed mass, G = 1";
  final parameter Real centralPotentialExact = -1.5 * GM / radius "-GM(3R^2 - r^2)/(2R^3) at r = 0";
  final parameter Real summary[4] = finiteDifferenceSummary(rho, radius, n, halfWidth);
  final parameter Real centralPotentialFD = summary[1] "Finite-difference potential at the centre node";
  final parameter Real sweeps = summary[2] "SOR sweeps used";
  final parameter Real residualFD = summary[3];
  final parameter Real convergedFD = summary[4] "1 when converged";
  final parameter Real rNode[nElements + 1] = {4.0 * radius * (i - 1) / nElements for i in 1:nElements + 1};
  final parameter Real rhoNode[nElements + 1] = {if rNode[i] < radius then rho elseif rNode[i] > radius then 0.0 else 0.5 * rho for i in 1:nElements + 1}
    "Half density on the surface node centres the linear ramp on the discontinuity";
  final parameter Real phiFE[nElements + 1] = GeneralRelativity.Numerical.poissonFiniteElementSpherical(rNode, rhoNode, 1.0);
  final parameter Real centralPotentialFE = phiFE[1];
  final parameter Real surfacePotentialFE = phiFE[integer(nElements / 4) + 1] "Node at r = radius";
  final parameter Real g[4, 4] = GeneralRelativity.Perturbation.weakFieldMetric(centralPotentialFE, zeros(3), c) "Weak-field metric at the centre";
  final parameter Real potentialRoundTrip = GeneralRelativity.Perturbation.newtonianPotential(g, c);
  Real finiteDifferenceError "Relative error of the finite-difference central potential";
  Real finiteElementError "Relative error of the finite-element central potential";
equation
  finiteDifferenceError = abs(centralPotentialFD / centralPotentialExact - 1.0);
  finiteElementError = abs(centralPotentialFE / centralPotentialExact - 1.0);
  assert(convergedFD > 0.5, "The finite-difference Poisson solve did not converge");
  assert(finiteDifferenceError < 0.05, "Finite-difference central potential is off by more than five percent");
  assert(finiteElementError < 5.0e-3, "Finite-element central potential is off by more than 0.5 percent");
  assert(abs(surfacePotentialFE / (-GM / radius) - 1.0) < 5.0e-3, "Finite-element surface potential disagrees with -GM/R");
  assert(abs(potentialRoundTrip - centralPotentialFE) < 1.0e-9 * abs(centralPotentialFE),
    "weakFieldMetric and newtonianPotential are not inverse to each other");
protected
  function finiteDifferenceSummary "Solve the box problem once and return {phi at centre, sweeps, residual, converged}"
    input Real rho;
    input Real radius;
    input Integer n;
    input Real halfWidth;
    output Real summary[4];
  protected
    GeneralRelativity.Numerical.Grid3D grid(n = {n, n, n}, origin = {-halfWidth, -halfWidth, -halfWidth},
      spacing = {2.0 * halfWidth / (n - 1), 2.0 * halfWidth / (n - 1), 2.0 * halfWidth / (n - 1)});
    Real GM = 4.0 / 3.0 * GeneralRelativity.Constants.pi * rho * radius^3;
    Real density[n, n, n];
    Real boundary[n, n, n];
    Real phi[n, n, n];
    Integer sweeps;
    Real residual;
    Boolean converged;
    Integer centre = integer((n + 1) / 2);
  algorithm
    density := GeneralRelativity.Numerical.uniformSphereDensity(grid, rho, radius, zeros(3));
    boundary := GeneralRelativity.Numerical.monopoleBoundary(grid, GM, zeros(3));
    (phi, sweeps, residual, converged) := GeneralRelativity.Numerical.poissonFiniteDifference3D(
      density, grid.spacing, boundary, 1.0, 2.0 / (1.0 + sin(GeneralRelativity.Constants.pi / n)), 4000, 1.0e-10);
    summary := {phi[centre, centre, centre], sweeps, residual, if converged then 1.0 else 0.0};
  end finiteDifferenceSummary;
  annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 1.0e-8, Interval = 0.1),
    Documentation(info="<html>
<p>Both solvers are given the same homogeneous sphere. The finite-difference
solution on a 61&sup3; box with the monopole far field as the boundary value is
limited by how the sphere is sampled on the grid (a few percent); the
spherical finite-element solution with one hundred linear elements is limited
by the density jump at the surface to a few parts in <code>10<sup>3</sup></code>;
a continuous density profile is nodally accurate to <code>10<sup>-6</sup></code>. The centre potential is lifted to the weak-field metric
and read back, closing the loop between
<code>&nabla;&sup2;&Phi; = 4&pi;G&rho;</code> and <code>G<sub>00</sub> = 8&pi;G T<sub>00</sub>/c<sup>4</sup></code>.</p>
</html>"));
end WeakFieldPoisson;
