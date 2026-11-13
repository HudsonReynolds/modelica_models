within GeneralRelativity.Connections;
function derivative "Central-difference connection gradient dGamma[m, i, j, k] = d Gamma^i_jk / d x^m"
  input GeneralRelativity.Interfaces.ConnectionFunction connection "Connection coefficients";
  input Real x[4] "Coordinates";
  input Real p[:] "Connection parameters";
  input Real step[4] "Finite-difference step per coordinate";
  output Real dGamma[4, 4, 4, 4];
protected
  Real xPlus[4];
  Real xMinus[4];
  Real GammaPlus[4, 4, 4];
  Real GammaMinus[4, 4, 4];
algorithm
  for m in 1:4 loop
    xPlus := x;
    xMinus := x;
    xPlus[m] := x[m] + step[m];
    xMinus[m] := x[m] - step[m];
    // Separate calls: OpenModelica's C backend cannot combine two functional-input
    // calls inside one slice assignment.
    GammaPlus := connection(xPlus, p);
    GammaMinus := connection(xMinus, p);
    // Explicit loops: a 4-D slice assignment is not reliable in OpenModelica's C backend.
    for i in 1:4 loop
      for j in 1:4 loop
        for k in 1:4 loop
          dGamma[m, i, j, k] := (GammaPlus[i, j, k] - GammaMinus[i, j, k]) / (2.0 * step[m]);
        end for;
      end for;
    end for;
  end for;
end derivative;
