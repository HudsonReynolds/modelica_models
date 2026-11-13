within GeneralRelativity.ExteriorCalculus;
function interior4 "Interior product of a vector with a 4-form, (i_X V)_bcd = X^a V_abcd"
  input Real X[4] "Vector X^a";
  input Real V[4, 4, 4, 4] "4-form";
  output Real H[4, 4, 4] "3-form";
algorithm
  for b in 1:4 loop
    for c in 1:4 loop
      for d in 1:4 loop
        H[b, c, d] := sum(X[a] * V[a, b, c, d] for a in 1:4);
      end for;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p><code>i<sub>X</sub>&epsilon; = *X<sup>&flat;</sup></code>: the volume form contracted
with a vector is the Hodge dual of the lowered vector.</p>
</html>"));
end interior4;
