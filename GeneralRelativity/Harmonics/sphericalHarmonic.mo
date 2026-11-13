within GeneralRelativity.Harmonics;
function sphericalHarmonic "Complex orthonormal spherical harmonic Y_lm(theta, phi) as real and imaginary parts"
  input Integer l "Degree, l >= 0";
  input Integer m "Order, -l <= m <= l";
  input Real theta(unit = "rad") "Polar angle from the z axis";
  input Real phi(unit = "rad") "Azimuth from the x axis";
  output Real re "Re Y_lm";
  output Real im "Im Y_lm";
protected
  Real P = normalizedLegendre(l, m, cos(theta));
algorithm
  re := P * cos(m * phi);
  im := P * sin(m * phi);
  annotation(Documentation(info="<html>
<p><code>Y<sub>lm</sub>(&theta;, &phi;) = sqrt((2l + 1)/(4&pi;) (l - m)!/(l + m)!) P<sub>l</sub><sup>m</sup>(cos&theta;) e<sup>im&phi;</sup></code>
with the Condon-Shortley phase, e.g.
<code>Y<sub>10</sub> = sqrt(3/(4&pi;)) cos&theta;</code> and
<code>Y<sub>11</sub> = -sqrt(3/(8&pi;)) sin&theta; e<sup>i&phi;</sup></code>. The addition theorem
<code>&Sigma;<sub>m</sub>|Y<sub>lm</sub>|&sup2; = (2l + 1)/(4&pi;)</code> holds at every point.</p>
</html>"));
end sphericalHarmonic;
