within GeneralRelativity.ExteriorCalculus;
function wedge22 "Wedge product of two 2-forms, a 4-form"
  input Real F[4, 4] "2-form (antisymmetric)";
  input Real G[4, 4] "2-form (antisymmetric)";
  output Real V[4, 4, 4, 4] "4-form F ^ G";
algorithm
  for i in 1:4 loop
    for j in 1:4 loop
      for k in 1:4 loop
        for l in 1:4 loop
          V[i, j, k, l] := F[i, j] * G[k, l] - F[i, k] * G[j, l] + F[i, l] * G[j, k]
                         + F[j, k] * G[i, l] - F[j, l] * G[i, k] + F[k, l] * G[i, j];
        end for;
      end for;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>The six (2,2) shuffles with their signs; symmetric under <code>F &harr; G</code>, and
<code>(a&and;b)&and;(c&and;d) = det[a; b; c; d] dx<sup>0</sup>&and;dx<sup>1</sup>&and;dx<sup>2</sup>&and;dx<sup>3</sup></code>.</p>
</html>"));
end wedge22;
