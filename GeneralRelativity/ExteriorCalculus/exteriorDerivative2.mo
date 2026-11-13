within GeneralRelativity.ExteriorCalculus;
function exteriorDerivative2 "Exterior derivative of a 2-form field, (dF)_abc = d_a F_bc + d_b F_ca + d_c F_ab"
  input TwoFormField field "2-form field";
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Field parameters";
  input Real step[4] "Finite-difference step per coordinate";
  output Real H[4, 4, 4] "3-form dF";
protected
  Real xPlus[4];
  Real xMinus[4];
  Real FPlus[4, 4];
  Real FMinus[4, 4];
  Real dPart[4, 4, 4] "dPart[a, b, c] = d_a F_bc";
algorithm
  for a in 1:4 loop
    xPlus := x;
    xMinus := x;
    xPlus[a] := x[a] + step[a];
    xMinus[a] := x[a] - step[a];
    FPlus := field(xPlus, p);
    FMinus := field(xMinus, p);
    dPart[a, :, :] := (FPlus - FMinus) / (2.0 * step[a]);
  end for;
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        H[a, b, c] := dPart[a, b, c] + dPart[b, c, a] + dPart[c, a, b];
      end for;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p><code>dF = 0</code> is the homogeneous Maxwell equation (no magnetic monopoles,
Faraday's law).</p>
</html>"));
end exteriorDerivative2;
