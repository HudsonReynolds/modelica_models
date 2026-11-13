within GeneralRelativity.ExteriorCalculus;
function interior2 "Interior product of a vector with a 2-form, (i_X F)_b = X^a F_ab"
  input Real X[4] "Vector X^a";
  input Real F[4, 4] "2-form";
  output Real w[4] "1-form";
algorithm
  w := X * F;
  annotation(Documentation(info="<html>
<p>Contraction on the first slot. With the wedge normalization of this package it
is an antiderivation: <code>i<sub>X</sub>(a&and;b) = (i<sub>X</sub>a) b - (i<sub>X</sub>b) a</code>.
For the Faraday tensor, <code>-i<sub>u</sub>F</code> is the electric field
measured by an observer with four-velocity <code>u</code> (up to the sign convention of <code>F</code>).</p>
</html>"));
end interior2;
