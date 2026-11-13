within GeneralRelativity.ExteriorCalculus;
function volumeElement "Metric volume density sqrt(|det g|) of a 4-metric"
  input Real g[4, 4] "Covariant metric";
  output Real sqrtDet "sqrt(|det g|)";
protected
  Real det;
  Integer s;
algorithm
  det := 0.0;
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        for d in 1:4 loop
          s := leviCivitaSymbol(a, b, c, d);
          if s <> 0 then
            det := det + s * g[1, a] * g[2, b] * g[3, c] * g[4, d];
          end if;
        end for;
      end for;
    end for;
  end for;
  sqrtDet := sqrt(abs(det));
  annotation(Documentation(info="<html>
<p>Leibniz expansion of the determinant (24 terms), no external linear algebra.
In a spherical chart of Schwarzschild <code>&radic;|g| = r&sup2; sin&theta;</code>.</p>
</html>"));
end volumeElement;
