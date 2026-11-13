within GeneralRelativity.Harmonics;
function realSphericalHarmonic "Real orthonormal spherical harmonic: cos(m phi) for m > 0, sin(|m| phi) for m < 0"
  input Integer l "Degree, l >= 0";
  input Integer m "Order, -l <= m <= l";
  input Real theta(unit = "rad") "Polar angle from the z axis";
  input Real phi(unit = "rad") "Azimuth from the x axis";
  output Real y "Real harmonic";
protected
  Integer mAbs = abs(m);
  Real P = normalizedLegendre(l, mAbs, cos(theta));
  Real sign = if mod(mAbs, 2) == 0 then 1.0 else -1.0 "(-1)^m, cancelling the Condon-Shortley phase";
algorithm
  if m > 0 then
    y := sqrt(2.0) * sign * P * cos(m * phi);
  elseif m < 0 then
    y := sqrt(2.0) * sign * P * sin(mAbs * phi);
  else
    y := P;
  end if;
  annotation(Documentation(info="<html>
<p><code>Y<sub>lm</sub><sup>real</sup> = &radic;2 (-1)<sup>m</sup> Re Y<sub>lm</sub></code> for
<code>m &gt; 0</code>, <code>Y<sub>l0</sub></code> for <code>m = 0</code> and
<code>&radic;2 (-1)<sup>m</sup> Im Y<sub>l|m|</sub></code> for <code>m &lt; 0</code>. The factor
<code>(-1)<sup>m</sup></code> removes the Condon-Shortley phase, so the real harmonics
are positive multiples of <code>cos(m&phi;)</code> and <code>sin(|m|&phi;)</code>:
<code>Y<sub>1,1</sub> = sqrt(3/(4&pi;)) x/r</code>, <code>Y<sub>1,-1</sub> = sqrt(3/(4&pi;)) y/r</code>.
They are orthonormal on the unit sphere and are the basis of geodesy gravity models
(up to the <code>4&pi;</code> normalization used there).</p>
</html>"));
end realSphericalHarmonic;
