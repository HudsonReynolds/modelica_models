within GeneralRelativity.Interfaces;
partial function PotentialFunction "Newtonian potential Phi(x) of a weak-field source"
  input Real x[3](each unit = "m") "Cartesian position";
  input Real p[:] "Potential parameter vector; layout is documented by each potential";
  output Real phi(unit = "m2/s2") "Newtonian potential, negative near mass";
end PotentialFunction;
