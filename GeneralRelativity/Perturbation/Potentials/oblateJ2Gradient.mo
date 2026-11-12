within GeneralRelativity.Perturbation.Potentials;
function oblateJ2Gradient "Gradient of oblateJ2; p = {GM, R, J2}"
  extends GeneralRelativity.Interfaces.PotentialGradientFunction;
protected
  Real GM = p[1];
  Real R = p[2];
  Real J2 = p[3];
  Real r2 = x * x;
  Real r = sqrt(r2);
  Real z2r2 = x[3]^2 / r2;
  Real k = 1.5 * J2 * GM * R^2 / (r2^2 * r);
algorithm
  // Newtonian point mass plus the standard J2 acceleration written as a gradient.
  gradient := GM * x / r^3;
  gradient[1] := gradient[1] + k * x[1] * (1.0 - 5.0 * z2r2);
  gradient[2] := gradient[2] + k * x[2] * (1.0 - 5.0 * z2r2);
  gradient[3] := gradient[3] + k * x[3] * (3.0 - 5.0 * z2r2);
  annotation(Documentation(info="<html>
<p>Minus this gradient is the familiar J2 perturbing acceleration
<code>-(3/2) J2 GM R&sup2;/r<sup>5</sup> [x(1 - 5z&sup2;/r&sup2;), y(1 - 5z&sup2;/r&sup2;), z(3 - 5z&sup2;/r&sup2;)]</code>
on top of <code>-GM x/r&sup3;</code>.</p>
</html>"));
end oblateJ2Gradient;
