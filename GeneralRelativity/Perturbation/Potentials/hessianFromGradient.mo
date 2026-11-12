within GeneralRelativity.Perturbation.Potentials;
function hessianFromGradient "Symmetrized central difference of an analytic potential gradient"
  input GeneralRelativity.Interfaces.PotentialGradientFunction gradient "Gradient function";
  input Real x[3](each unit = "m") "Position";
  input Real p[:] "Potential parameters";
  input Real step(unit = "m") "Finite-difference step";
  output Real hessian[3, 3](each unit = "1/s2");
protected
  Real xPlus[3];
  Real xMinus[3];
  Real gPlus[3];
  Real gMinus[3];
  Real column[3];
algorithm
  for j in 1:3 loop
    xPlus := x;
    xMinus := x;
    xPlus[j] := x[j] + step;
    xMinus[j] := x[j] - step;
    gPlus := gradient(xPlus, p);
    gMinus := gradient(xMinus, p);
    column := (gPlus - gMinus) / (2.0 * step);
    hessian[:, j] := column;
  end for;
  hessian := 0.5 * (hessian + transpose(hessian));
end hessianFromGradient;
