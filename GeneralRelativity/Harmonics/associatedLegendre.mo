within GeneralRelativity.Harmonics;
function associatedLegendre "Associated Legendre function P_l^m(x) with the Condon-Shortley phase"
  input Integer l "Degree, l >= 0";
  input Integer m "Order, -l <= m <= l (zero is returned otherwise)";
  input Real x "Argument in [-1, 1]";
  output Real P "P_l^m(x)";
protected
  Integer mAbs = abs(m);
  Real s = sqrt(max(1.0 - x^2, 0.0));
  Real pmm "P_m^m";
  Real pPrev;
  Real pCur;
  Real pNext;
algorithm
  if l < 0 or mAbs > l then
    P := 0.0;
  else
    // P_m^m = (-1)^m (2m - 1)!! (1 - x^2)^(m/2)
    pmm := 1.0;
    for k in 1:mAbs loop
      pmm := -pmm * (2 * k - 1) * s;
    end for;
    if l == mAbs then
      P := pmm;
    else
      pPrev := pmm;
      pCur := x * (2 * mAbs + 1) * pmm;
      for ll in mAbs + 2:l loop
        pNext := (x * (2 * ll - 1) * pCur - (ll + mAbs - 1) * pPrev) / (ll - mAbs);
        pPrev := pCur;
        pCur := pNext;
      end for;
      P := pCur;
    end if;
    if m < 0 then
      P := (if mod(mAbs, 2) == 0 then 1.0 else -1.0) * exp(logFactorial(l - mAbs) - logFactorial(l + mAbs)) * P;
    end if;
  end if;
  annotation(Documentation(info="<html>
<p><code>P<sub>l</sub><sup>m</sup>(x) = (-1)<sup>m</sup>(1 - x&sup2;)<sup>m/2</sup> d<sup>m</sup>P<sub>l</sub>(x)/dx<sup>m</sup></code>
by upward recurrence in <code>l</code> from
<code>P<sub>m</sub><sup>m</sup> = (-1)<sup>m</sup>(2m - 1)!!(1 - x&sup2;)<sup>m/2</sup></code>:
<code>(l - m)P<sub>l</sub><sup>m</sup> = x(2l - 1)P<sub>l-1</sub><sup>m</sup> - (l + m - 1)P<sub>l-2</sub><sup>m</sup></code>.
Negative orders use
<code>P<sub>l</sub><sup>-m</sup> = (-1)<sup>m</sup>(l - m)!/(l + m)! P<sub>l</sub><sup>m</sup></code>.
Examples: <code>P<sub>1</sub><sup>1</sup> = -sqrt(1 - x&sup2;)</code>,
<code>P<sub>2</sub><sup>1</sup> = -3x sqrt(1 - x&sup2;)</code>,
<code>P<sub>2</sub><sup>2</sup> = 3(1 - x&sup2;)</code>. The values grow like
<code>(2l - 1)!!</code>; for high degree use
<a href=\"modelica://GeneralRelativity.Harmonics.normalizedLegendre\">normalizedLegendre</a>.</p>
</html>"));
end associatedLegendre;
