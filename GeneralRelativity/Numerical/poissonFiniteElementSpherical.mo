within GeneralRelativity.Numerical;
function poissonFiniteElementSpherical "Linear finite elements for (r^2 phi')' = 4 pi G r^2 rho with zero flux at the centre and the monopole value at the outer node"
  input Real r[:](each unit = "m") "Node radii, increasing, r[1] >= 0";
  input Real rho[size(r, 1)](each unit = "kg/m3") "Mass density at the nodes, interpolated linearly";
  input Real gravitationalConstant(unit = "m3/(kg.s2)") = 1.0 "Newton's constant";
  output Real phi[size(r, 1)](each unit = "m2/s2") "Potential at the nodes";
  output Real enclosedMass[size(r, 1)](each unit = "kg") "Mass inside each node radius";
protected
  Integer n = size(r, 1);
  constant Real gaussPoint[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
  constant Real gaussWeight[3] = {0.5555555555555556, 0.8888888888888888, 0.5555555555555556};
  Real fourPiG = 4.0 * GeneralRelativity.Constants.pi * gravitationalConstant;
  Real lower[size(r, 1)] "Sub-diagonal";
  Real diag[size(r, 1)];
  Real upper[size(r, 1)] "Super-diagonal";
  Real load[size(r, 1)];
  Real length;
  Real rq;
  Real rhoq;
  Real na;
  Real nb;
  Real stiffness;
  Real shellMass;
  Real factor;
algorithm
  assert(n >= 2, "poissonFiniteElementSpherical needs at least two nodes");
  diag := zeros(n);
  lower := zeros(n);
  upper := zeros(n);
  load := zeros(n);
  enclosedMass[1] := 0.0;
  for e in 1:n - 1 loop
    length := r[e + 1] - r[e];
    assert(length > 0.0, "poissonFiniteElementSpherical needs strictly increasing radii");
    // Element stiffness: (1/L^2) * int r^2 dr over the element.
    stiffness := (r[e + 1]^3 - r[e]^3) / (3.0 * length^2);
    diag[e] := diag[e] + stiffness;
    diag[e + 1] := diag[e + 1] + stiffness;
    upper[e] := upper[e] - stiffness;
    lower[e + 1] := lower[e + 1] - stiffness;
    // Load -4 pi G int r^2 rho N_i dr by three-point Gauss quadrature, exact for this degree.
    shellMass := 0.0;
    for q in 1:3 loop
      nb := 0.5 * (1.0 + gaussPoint[q]);
      na := 1.0 - nb;
      rq := na * r[e] + nb * r[e + 1];
      rhoq := na * rho[e] + nb * rho[e + 1];
      factor := 0.5 * length * gaussWeight[q] * rq^2 * rhoq;
      load[e] := load[e] - fourPiG * factor * na;
      load[e + 1] := load[e + 1] - fourPiG * factor * nb;
      shellMass := shellMass + 4.0 * GeneralRelativity.Constants.pi * factor;
    end for;
    enclosedMass[e + 1] := enclosedMass[e] + shellMass;
  end for;
  // Dirichlet outer node: phi(R) = -G M / R.
  diag[n] := 1.0;
  lower[n] := 0.0;
  load[n] := -gravitationalConstant * enclosedMass[n] / r[n];
  // Thomas algorithm.
  for i in 2:n loop
    factor := lower[i] / diag[i - 1];
    diag[i] := diag[i] - factor * upper[i - 1];
    load[i] := load[i] - factor * load[i - 1];
  end for;
  phi[n] := load[n] / diag[n];
  for back in 1:n - 1 loop
    phi[n - back] := (load[n - back] - upper[n - back] * phi[n - back + 1]) / diag[n - back];
  end for;
  annotation(Documentation(info="<html>
<p>Galerkin weak form with hat functions. The stiffness matrix is exact for
<code>r&sup2;</code> weights and the load uses three-point Gauss quadrature, which is
exact for the quartic integrand. The zero-flux condition at <code>r = 0</code> is
natural, so no equation is modified there. Accuracy is second order in the
element length; a homogeneous sphere with fifty elements reproduces the central
potential <code>-3GM/(2R)</code> to about <code>10<sup>-4</sup></code>.</p>
</html>"));
end poissonFiniteElementSpherical;
