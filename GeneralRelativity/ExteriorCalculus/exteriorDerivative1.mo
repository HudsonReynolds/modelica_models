within GeneralRelativity.ExteriorCalculus;
function exteriorDerivative1 "Exterior derivative of a 1-form field, (dw)_ab = d_a w_b - d_b w_a"
  input OneFormField field "1-form field";
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Field parameters";
  input Real step[4] "Finite-difference step per coordinate";
  output Real F[4, 4] "2-form dw";
protected
  Real xPlus[4];
  Real xMinus[4];
  Real wPlus[4];
  Real wMinus[4];
  Real dPart[4, 4] "dPart[a, b] = d_a w_b";
algorithm
  for a in 1:4 loop
    xPlus := x;
    xMinus := x;
    xPlus[a] := x[a] + step[a];
    xMinus[a] := x[a] - step[a];
    wPlus := field(xPlus, p);
    wMinus := field(xMinus, p);
    dPart[a, :] := (wPlus - wMinus) / (2.0 * step[a]);
  end for;
  F := dPart - transpose(dPart);
  annotation(Documentation(info="<html>
<p>Partial derivatives suffice: the exterior derivative needs no connection. For
the electromagnetic potential, <code>F = dA</code> is the Faraday tensor.</p>
</html>"));
end exteriorDerivative1;
