within GeneralRelativity.Perturbation.Potentials;
function pointMassHessian "Phi_,ij = GM (delta_ij r^2 - 3 x_i x_j) / r^5, the traceless point-mass tidal tensor"
  input Real x[3](each unit = "m") "Position";
  input Real GM(unit = "m3/s2") "Gravitational parameter";
  output Real hessian[3, 3](each unit = "1/s2");
protected
  Real r2 = x * x;
  Real r5 = r2^2 * sqrt(r2);
algorithm
  for i in 1:3 loop
    for j in 1:3 loop
      hessian[i, j] := GM * ((if i == j then r2 else 0.0) - 3.0 * x[i] * x[j]) / r5;
    end for;
  end for;
end pointMassHessian;
