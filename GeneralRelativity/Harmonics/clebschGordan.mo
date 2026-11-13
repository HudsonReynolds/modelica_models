within GeneralRelativity.Harmonics;
function clebschGordan "Clebsch-Gordan coefficient <j1 m1 j2 m2 | J M> by the Racah formula; integer or half-integer arguments"
  input Real j1 "First angular momentum (integer or half-integer)";
  input Real m1 "Projection of j1";
  input Real j2 "Second angular momentum";
  input Real m2 "Projection of j2";
  input Real J "Total angular momentum";
  input Real M "Projection of J";
  output Real C "<j1 m1 j2 m2 | J M>; zero when a selection rule fails";
protected
  Integer tj1 = integer(floor(2.0 * j1 + 0.5)) "2 j1";
  Integer tm1 = integer(floor(2.0 * m1 + 0.5));
  Integer tj2 = integer(floor(2.0 * j2 + 0.5));
  Integer tm2 = integer(floor(2.0 * m2 + 0.5));
  Integer tJ = integer(floor(2.0 * J + 0.5));
  Integer tM = integer(floor(2.0 * M + 0.5));
  Integer a1;
  Integer a2;
  Integer a3;
  Integer a4;
  Integer b1;
  Integer b2;
  Integer kMin;
  Integer kMax;
  Real logPrefactor;
  Real total;
algorithm
  C := 0.0;
  if tj1 >= 0 and tj2 >= 0 and tJ >= 0 and tm1 + tm2 == tM
    and abs(tm1) <= tj1 and abs(tm2) <= tj2 and abs(tM) <= tJ
    and mod(tj1 + tm1, 2) == 0 and mod(tj2 + tm2, 2) == 0 and mod(tJ + tM, 2) == 0
    and tJ >= abs(tj1 - tj2) and tJ <= tj1 + tj2 and mod(tj1 + tj2 + tJ, 2) == 0 then
    a1 := div(tJ + tj1 - tj2, 2); // J + j1 - j2
    a2 := div(tJ - tj1 + tj2, 2); // J - j1 + j2
    a3 := div(tj1 + tj2 - tJ, 2); // j1 + j2 - J
    a4 := div(tj1 + tj2 + tJ, 2) + 1; // j1 + j2 + J + 1
    b1 := div(tJ - tj2 + tm1, 2); // J - j2 + m1
    b2 := div(tJ - tj1 - tm2, 2); // J - j1 - m2
    logPrefactor := 0.5 * (log(tJ + 1.0) + logFactorial(a1) + logFactorial(a2) + logFactorial(a3) - logFactorial(a4)
      + logFactorial(div(tJ + tM, 2)) + logFactorial(div(tJ - tM, 2))
      + logFactorial(div(tj1 - tm1, 2)) + logFactorial(div(tj1 + tm1, 2))
      + logFactorial(div(tj2 - tm2, 2)) + logFactorial(div(tj2 + tm2, 2)));
    kMin := max(0, max(-b1, -b2));
    kMax := min(a3, min(div(tj1 - tm1, 2), div(tj2 + tm2, 2)));
    total := 0.0;
    for k in kMin:kMax loop
      total := total + (if mod(k, 2) == 0 then 1.0 else -1.0) * exp(logPrefactor - logFactorial(k) - logFactorial(a3 - k)
        - logFactorial(div(tj1 - tm1, 2) - k) - logFactorial(div(tj2 + tm2, 2) - k)
        - logFactorial(b1 + k) - logFactorial(b2 + k));
    end for;
    C := total;
  end if;
  annotation(Documentation(info="<html>
<p>Racah's closed form</p>
<pre>&lang;j1 m1 j2 m2|J M&rang; = &delta;(M, m1 + m2)
  sqrt((2J+1)(J+j1-j2)!(J-j1+j2)!(j1+j2-J)!/(j1+j2+J+1)!)
  sqrt((J+M)!(J-M)!(j1-m1)!(j1+m1)!(j2-m2)!(j2+m2)!)
  &Sigma;<sub>k</sub> (-1)<sup>k</sup>/[k!(j1+j2-J-k)!(j1-m1-k)!(j2+m2-k)!(J-j2+m1+k)!(J-j1-m2+k)!]</pre>
<p>with the factorials evaluated as
<a href=\"modelica://GeneralRelativity.Harmonics.logFactorial\">logFactorial</a> and each term
exponentiated after the prefactor is folded in, so nothing overflows. Arguments are
<code>Real</code>; each is doubled and rounded, which makes half-integers such as
<code>0.5</code> exact. The coefficient is zero unless <code>M = m1 + m2</code>,
<code>|m| &le; j</code>, <code>j &plusmn; m</code> are integers and
<code>|j1 - j2| &le; J &le; j1 + j2</code> with <code>j1 + j2 + J</code> an integer. Example:
<code>&lang;&frac12; &frac12; &frac12; -&frac12;|0 0&rang; = 1/&radic;2</code>.</p>
</html>"));
end clebschGordan;
