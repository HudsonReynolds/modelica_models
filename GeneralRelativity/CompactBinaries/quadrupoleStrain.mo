within GeneralRelativity.CompactBinaries;
function quadrupoleStrain "Quadrupole gravitational-wave strain of a binary for an observer on the orbital axis, scaled by distance over mass"
  input Real x[3] "Relative position x1 - x2";
  input Real v[3] "Relative velocity";
  input Real M = 1.0 "Total mass";
  input Real eta = 0.25 "Symmetric mass ratio";
  output Real h[2] "{h+, hx} times R/M (R the observer's distance)";
protected
  Real r = sqrt(x * x);
algorithm
  h[1] := 2.0 * eta * (v[1]^2 - v[2]^2 - M * (x[1]^2 - x[2]^2) / r^3);
  h[2] := 4.0 * eta * (v[1] * v[2] - M * x[1] * x[2] / r^3);
  annotation(Documentation(info="<html>
<p>The quadrupole formula <code>h<sub>ij</sub> = (2/R) d&sup2;I<sub>ij</sub>/dt&sup2;</code> with
<code>I<sub>ij</sub> = &mu; x<sub>i</sub>x<sub>j</sub></code> and the Newtonian acceleration
<code>&minus;Mx/r&sup3;</code>, for an observer on the <code>+z</code> axis (face-on):
<code>h<sub>+</sub> = (h<sub>xx</sub> &minus; h<sub>yy</sub>)/2</code>, <code>h<sub>&times;</sub> = h<sub>xy</sub></code>.
On a circular orbit the amplitude is <code>R|h|/M = 4&eta;(M&Omega;)<sup>2/3</sup></code> and the wave
frequency is twice the orbital one. Multiply by <code>M/R</code> (both in metres) for the strain.</p>
</html>"));
end quadrupoleStrain;
