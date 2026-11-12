within GeneralRelativity.Perturbation;
function weakFieldCurvature "Linearized Riemann and Einstein tensors of a static Newtonian potential from its Hessian"
  input Real hessian[3, 3] "d^2 Phi / (d x_i d x_j)";
  input Real c(unit = "m/s") = 1.0 "Speed of light";
  output Real R[4, 4, 4, 4] "Linearized R_abcd";
  output Real G[4, 4] "Linearized G_ab";
  output Real laplacian(unit = "1/s2") "trace of the Hessian, 4 pi G rho for a Newtonian source";
protected
  Real ddh[4, 4, 4, 4];
  Real Ric[4, 4];
  Real scalar;
algorithm
  ddh := zeros(4, 4, 4, 4);
  // h_00 = h_ii = -2 Phi / c^2, static, so only spatial second derivatives survive.
  for i in 1:3 loop
    for j in 1:3 loop
      ddh[i + 1, j + 1, 1, 1] := -2.0 * hessian[i, j] / c^2;
      for k in 1:3 loop
        ddh[i + 1, j + 1, k + 1, k + 1] := -2.0 * hessian[i, j] / c^2;
      end for;
    end for;
  end for;
  R := linearizedRiemann(ddh);
  (G, Ric, scalar) := linearizedEinstein(ddh);
  laplacian := hessian[1, 1] + hessian[2, 2] + hessian[3, 3];
end weakFieldCurvature;
