within GeneralRelativity.Numerical;
function poissonFiniteDifference3D "Solve laplacian(phi) = 4 pi G rho on a Cartesian grid by seven-point successive over-relaxation with Dirichlet boundaries"
  input Real rho[:, :, :](each unit = "kg/m3") "Mass density at the nodes";
  input Real spacing[3](each unit = "m") "Node spacing per axis";
  input Real phiBoundary[size(rho, 1), size(rho, 2), size(rho, 3)](each unit = "m2/s2")
    "Dirichlet values; only boundary nodes are read, interior entries seed the iteration";
  input Real gravitationalConstant(unit = "m3/(kg.s2)") = 1.0 "Newton's constant";
  input Real relaxation = 1.5 "SOR factor in (0, 2); 1 is Gauss-Seidel";
  input Integer maxIterations = 2000 "Maximum number of full sweeps";
  input Real tolerance(unit = "m2/s2") = 1.0e-9 "Stop when the largest nodal update is below this";
  output Real phi[size(rho, 1), size(rho, 2), size(rho, 3)](each unit = "m2/s2") "Potential";
  output Integer iterations "Sweeps performed";
  output Real residual(unit = "m2/s2") "Largest nodal update of the last sweep";
  output Boolean converged;
protected
  Integer nx = size(rho, 1);
  Integer ny = size(rho, 2);
  Integer nz = size(rho, 3);
  Real ax = 1.0 / spacing[1]^2;
  Real ay = 1.0 / spacing[2]^2;
  Real az = 1.0 / spacing[3]^2;
  Real diagonal = 2.0 * (ax + ay + az);
  Real source = 4.0 * GeneralRelativity.Constants.pi * gravitationalConstant;
  Real gaussSeidel;
  Real update;
algorithm
  assert(nx >= 3 and ny >= 3 and nz >= 3, "poissonFiniteDifference3D needs at least three nodes per axis");
  phi := phiBoundary;
  iterations := 0;
  residual := 0.0;
  converged := false;
  for sweep in 1:maxIterations loop
    if not converged then
      residual := 0.0;
      for i in 2:nx - 1 loop
        for j in 2:ny - 1 loop
          for k in 2:nz - 1 loop
            gaussSeidel := (ax * (phi[i + 1, j, k] + phi[i - 1, j, k])
              + ay * (phi[i, j + 1, k] + phi[i, j - 1, k])
              + az * (phi[i, j, k + 1] + phi[i, j, k - 1])
              - source * rho[i, j, k]) / diagonal;
            update := relaxation * (gaussSeidel - phi[i, j, k]);
            phi[i, j, k] := phi[i, j, k] + update;
            residual := max(residual, abs(update));
          end for;
        end for;
      end for;
      iterations := sweep;
      converged := residual <= tolerance;
    end if;
  end for;
  annotation(Documentation(info="<html>
<p>Second-order accurate in the spacing. The optimal SOR factor for an
<code>n</code>-node cube is near <code>2/(1 + sin(&pi;/n))</code>. With the monopole
far field as the boundary value the solution outside the source converges to
<code>-GM/r</code> and inside a homogeneous sphere to
<code>-GM(3R&sup2; - r&sup2;)/(2R&sup3;)</code>. The weak-field metric follows from
<a href=\"modelica://GeneralRelativity.Perturbation.weakFieldMetric\">Perturbation.weakFieldMetric</a>
at every node.</p>
</html>"));
end poissonFiniteDifference3D;
