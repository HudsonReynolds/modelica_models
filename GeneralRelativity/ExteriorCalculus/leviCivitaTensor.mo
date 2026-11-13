within GeneralRelativity.ExteriorCalculus;
function leviCivitaTensor "Covariant Levi-Civita tensor (volume form) eps_abcd = sqrt(|det g|) [abcd]"
  input Real g[4, 4] "Covariant metric";
  output Real eps[4, 4, 4, 4] "eps_abcd, with eps_1234 = +sqrt(|det g|)";
protected
  Real root;
algorithm
  root := volumeElement(g);
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        for d in 1:4 loop
          eps[a, b, c, d] := root * leviCivitaSymbol(a, b, c, d);
        end for;
      end for;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>Orientation: <code>&epsilon;<sub>0123</sub> = +&radic;|g|</code> in physics labels
(<code>eps[1,2,3,4]</code>), so <code>&epsilon; = &radic;|g| dx<sup>0</sup>&and;dx<sup>1</sup>&and;dx<sup>2</sup>&and;dx<sup>3</sup></code>.
In Lorentzian signature the contravariant tensor has
<code>&epsilon;<sup>0123</sup> = -1/&radic;|g|</code> and
<code>&epsilon;<sub>abcd</sub>&epsilon;<sup>abcd</sup> = -24</code>.</p>
</html>"));
end leviCivitaTensor;
