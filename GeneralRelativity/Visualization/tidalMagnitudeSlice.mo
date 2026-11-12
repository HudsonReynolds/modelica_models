within GeneralRelativity.Visualization;
function tidalMagnitudeSlice "Frobenius norm of the Newtonian tidal tensor (Hessian of the potential) on a z = offset plane"
  input GeneralRelativity.Interfaces.PotentialGradientFunction gradient "Potential gradient";
  input Real p[:] "Potential parameters";
  input Real step(unit = "m") "Finite-difference step for the Hessian";
  input Real u[:](each unit = "m") "x samples";
  input Real v[:](each unit = "m") "y samples";
  input Real offset(unit = "m") = 0.0 "z of the plane";
  input Real innerRadius(unit = "m") = 0.0 "Points closer to the origin than this are evaluated on the sphere of this radius, which keeps a point-mass singularity out of the colour scale";
  output Real values[size(u, 1), size(v, 1)](each unit = "1/s2");
protected
  Real hessian[3, 3];
  Real point[3];
  Real distance;
algorithm
  for i in 1:size(u, 1) loop
    for j in 1:size(v, 1) loop
      point := {u[i], v[j], offset};
      distance := sqrt(point * point);
      if distance < innerRadius then
        point := if distance > 1.0e-9 * innerRadius then point * (innerRadius / distance) else {innerRadius, 0.0, 0.0};
      end if;
      hessian := GeneralRelativity.Perturbation.Potentials.hessianFromGradient(gradient, point, p, step);
      values[i, j] := sqrt(sum(hessian[a, b]^2 for a in 1:3, b in 1:3));
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>In the weak field the Riemann tensor is the tidal tensor
<code>R<sub>0i0j</sub> = &Phi;<sub>,ij</sub>/c&sup2;</code>, so this norm is
<code>c&sup2;</code> times the size of the curvature; it scales as
<code>GM/r&sup3;</code> for a point mass.</p>
</html>"));
end tidalMagnitudeSlice;
