within GeneralRelativity.Interfaces;
partial function PotentialGradientFunction "Gradient of a Newtonian potential; minus the gradient is the acceleration"
  input Real x[3](each unit = "m") "Cartesian position";
  input Real p[:] "Potential parameter vector; layout is documented by each potential";
  output Real gradient[3](each unit = "m/s2") "d Phi / d x";
end PotentialGradientFunction;
