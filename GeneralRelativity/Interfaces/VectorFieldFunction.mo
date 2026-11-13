within GeneralRelativity.Interfaces;
partial function VectorFieldFunction "Contravariant vector field xi^i(x), for example a candidate Killing vector"
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Field parameter vector; layout is documented by each field";
  output Real xi[4] "Contravariant components xi^i";
end VectorFieldFunction;
