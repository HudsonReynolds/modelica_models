within GeneralRelativity.ExteriorCalculus;
function wedge11 "Wedge product of two 1-forms, (a ^ b)_ij = a_i b_j - a_j b_i"
  input Real a[4] "1-form";
  input Real b[4] "1-form";
  output Real F[4, 4] "2-form a ^ b";
algorithm
  for i in 1:4 loop
    for j in 1:4 loop
      F[i, j] := a[i] * b[j] - a[j] * b[i];
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>Forms are stored as fully antisymmetric component arrays with
<code>&omega; = (1/p!) &omega;<sub>a<sub>1</sub>..a<sub>p</sub></sub> dx<sup>a<sub>1</sub></sup>&and;..&and;dx<sup>a<sub>p</sub></sup></code>
and <code>(&alpha;&and;&beta;)<sub>a..b..</sub> = ((p+q)!/(p!q!)) &alpha;<sub>[a..</sub>&beta;<sub>b..]</sub></code>;
see <a href=\"modelica://GeneralRelativity.ExteriorCalculus\">ExteriorCalculus</a>.</p>
</html>"));
end wedge11;
