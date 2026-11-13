within GeneralRelativity.Symmetry;
function fieldDerivative "Central-difference gradient dxi[a, c] = d xi^c / d x^a of a vector field"
  input GeneralRelativity.Interfaces.VectorFieldFunction field "Vector field";
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real q[:] "Field parameters";
  input Real step[4] "Finite-difference step per coordinate";
  output Real dxi[4, 4] "dxi[a, c] = d_a xi^c";
protected
  Real xPlus[4];
  Real xMinus[4];
  Real xiPlus[4];
  Real xiMinus[4];
algorithm
  for a in 1:4 loop
    xPlus := x;
    xMinus := x;
    xPlus[a] := x[a] + step[a];
    xMinus[a] := x[a] - step[a];
    xiPlus := field(xPlus, q);
    xiMinus := field(xMinus, q);
    dxi[a, :] := (xiPlus - xiMinus) / (2.0 * step[a]);
  end for;
end fieldDerivative;
