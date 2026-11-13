within GeneralRelativity.ExteriorCalculus;
function hodge1 "Hodge dual of a 1-form, (*w)_bcd = w^a eps_abcd"
  input Real w[4] "1-form w_a";
  input Real g[4, 4] "Covariant metric";
  output Real H[4, 4, 4] "3-form *w";
protected
  Real gInv[4, 4];
  Boolean accepted;
  Real eps[4, 4, 4, 4];
  Real wUp[4];
algorithm
  (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(g);
  eps := leviCivitaTensor(g);
  wUp := gInv * w;
  for b in 1:4 loop
    for c in 1:4 loop
      for d in 1:4 loop
        H[b, c, d] := sum(wUp[a] * eps[a, b, c, d] for a in 1:4);
      end for;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>In Minkowski space <code>*dx<sup>0</sup> = -dx<sup>1</sup>&and;dx<sup>2</sup>&and;dx<sup>3</sup></code>.</p>
</html>"));
end hodge1;
