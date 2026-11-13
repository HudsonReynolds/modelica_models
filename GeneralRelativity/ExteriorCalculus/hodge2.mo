within GeneralRelativity.ExteriorCalculus;
function hodge2 "Hodge dual of a 2-form, (*F)_cd = (1/2) F^ab eps_abcd"
  input Real F[4, 4] "2-form F_ab (antisymmetric)";
  input Real g[4, 4] "Covariant metric";
  output Real dual[4, 4] "2-form *F";
protected
  Real gInv[4, 4];
  Boolean accepted;
  Real eps[4, 4, 4, 4];
  Real FUp[4, 4];
algorithm
  (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(g);
  eps := leviCivitaTensor(g);
  FUp := gInv * F * gInv;
  for c in 1:4 loop
    for d in 1:4 loop
      dual[c, d] := 0.0;
      for a in 1:4 loop
        for b in 1:4 loop
          dual[c, d] := dual[c, d] + 0.5 * FUp[a, b] * eps[a, b, c, d];
        end for;
      end for;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>In Minkowski space <code>*(dx<sup>0</sup>&and;dx<sup>1</sup>) = -dx<sup>2</sup>&and;dx<sup>3</sup></code>
and <code>**F = -F</code>: the electromagnetic duality <code>(E, B) &rarr; (-B, E)</code>
squares to minus one.</p>
</html>"));
end hodge2;
