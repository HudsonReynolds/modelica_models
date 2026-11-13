within GeneralRelativity.Harmonics;
function normalizedLegendre "Orthonormalized associated Legendre function, Y_lm(theta, 0) as a function of x = cos(theta)"
  input Integer l "Degree, l >= 0";
  input Integer m "Order, -l <= m <= l (zero is returned otherwise)";
  input Real x "Argument cos(theta) in [-1, 1]";
  output Real P "sqrt((2l + 1)/(4 pi) (l - m)!/(l + m)!) P_l^m(x)";
protected
  Integer mAbs = abs(m);
  Real s = sqrt(max(1.0 - x^2, 0.0));
  Real pPrev;
  Real pCur;
  Real pNext;
algorithm
  if l < 0 or mAbs > l then
    P := 0.0;
  else
    // Sectoral seed: Pbar_m^m = -sqrt((2m + 1)/(2m)) s Pbar_(m-1)^(m-1), Pbar_0^0 = 1/sqrt(4 pi).
    pCur := 1.0 / sqrt(4.0 * GeneralRelativity.Constants.pi);
    for k in 1:mAbs loop
      pCur := -sqrt((2.0 * k + 1.0) / (2.0 * k)) * s * pCur;
    end for;
    if l > mAbs then
      pPrev := pCur;
      pCur := x * sqrt(2.0 * mAbs + 3.0) * pPrev;
      for ll in mAbs + 2:l loop
        pNext := sqrt((4.0 * ll^2 - 1.0) / (ll^2 - mAbs^2))
          * (x * pCur - sqrt(((ll - 1.0)^2 - mAbs^2) / (4.0 * (ll - 1.0)^2 - 1.0)) * pPrev);
        pPrev := pCur;
        pCur := pNext;
      end for;
    end if;
    P := if m < 0 and mod(mAbs, 2) == 1 then -pCur else pCur;
  end if;
  annotation(Documentation(info="<html>
<p>The normalization makes <code>Y<sub>lm</sub> = P e<sup>im&phi;</sup></code> orthonormal on the
unit sphere; the Condon-Shortley phase is included and
<code>P(l, -m) = (-1)<sup>m</sup>P(l, m)</code>. The three-term recurrence on the
normalized functions has coefficients of order one, so it neither overflows nor loses
accuracy for large <code>l</code>; only <code>sin<sup>m</sup>&theta;</code> can underflow very
near the poles at very high order.</p>
</html>"));
end normalizedLegendre;
