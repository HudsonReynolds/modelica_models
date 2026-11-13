within GeneralRelativity.ExteriorCalculus;
function exteriorDerivative3 "Exterior derivative of a 3-form field, a 4-form"
  input ThreeFormField field "3-form field";
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Field parameters";
  input Real step[4] "Finite-difference step per coordinate";
  output Real V[4, 4, 4, 4] "(dH)_abcd = d_a H_bcd - d_b H_acd + d_c H_abd - d_d H_abc";
protected
  Real xPlus[4];
  Real xMinus[4];
  Real HPlus[4, 4, 4];
  Real HMinus[4, 4, 4];
  Real dPart[4, 4, 4, 4] "dPart[a, b, c, d] = d_a H_bcd";
algorithm
  for a in 1:4 loop
    xPlus := x;
    xMinus := x;
    xPlus[a] := x[a] + step[a];
    xMinus[a] := x[a] - step[a];
    HPlus := field(xPlus, p);
    HMinus := field(xMinus, p);
    for b in 1:4 loop
      for c in 1:4 loop
        for d in 1:4 loop
          dPart[a, b, c, d] := (HPlus[b, c, d] - HMinus[b, c, d]) / (2.0 * step[a]);
        end for;
      end for;
    end for;
  end for;
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        for d in 1:4 loop
          V[a, b, c, d] := dPart[a, b, c, d] - dPart[b, a, c, d] + dPart[c, a, b, d] - dPart[d, a, b, c];
        end for;
      end for;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>For <code>H = *w</code>, <code>d*w = (&nabla;<sub>a</sub>w<sup>a</sup>) &epsilon;</code>.</p>
</html>"));
end exteriorDerivative3;
