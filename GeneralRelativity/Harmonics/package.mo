within GeneralRelativity;
package Harmonics "Associated Legendre functions, spherical harmonics, Clebsch-Gordan and Wigner 3j coefficients, Gaunt integrals"
  annotation(Documentation(info="<html>
<h4>Conventions</h4>
<ul>
<li><a href=\"modelica://GeneralRelativity.Harmonics.associatedLegendre\">associatedLegendre</a>
includes the Condon-Shortley phase:
<code>P<sub>l</sub><sup>m</sup>(x) = (-1)<sup>m</sup>(1 - x&sup2;)<sup>m/2</sup> d<sup>m</sup>P<sub>l</sub>/dx<sup>m</sup></code>,
so <code>P<sub>1</sub><sup>1</sup>(x) = -sqrt(1 - x&sup2;)</code>, and
<code>P<sub>l</sub><sup>-m</sup> = (-1)<sup>m</sup>(l - m)!/(l + m)! P<sub>l</sub><sup>m</sup></code>.</li>
<li><a href=\"modelica://GeneralRelativity.Harmonics.sphericalHarmonic\">sphericalHarmonic</a>
<code>Y<sub>lm</sub>(&theta;, &phi;) = sqrt((2l + 1)/(4&pi;) (l - m)!/(l + m)!) P<sub>l</sub><sup>m</sup>(cos&theta;) e<sup>im&phi;</sup></code>
is orthonormal on the unit sphere, <code>&int;Y<sub>lm</sub>Y*<sub>l'm'</sub>d&Omega; = &delta;<sub>ll'</sub>&delta;<sub>mm'</sub></code>,
and <code>Y<sub>l,-m</sub> = (-1)<sup>m</sup>Y*<sub>lm</sub></code> (Jackson, Sakurai, Arfken).</li>
<li><a href=\"modelica://GeneralRelativity.Harmonics.realSphericalHarmonic\">realSphericalHarmonic</a>
is the real orthonormal basis <code>&radic;2 (-1)<sup>m</sup> Re Y<sub>lm</sub></code> for
<code>m &gt; 0</code> and <code>&radic;2 (-1)<sup>m</sup> Im Y<sub>l|m|</sub></code> for <code>m &lt; 0</code>,
so <code>Y<sub>1,1</sub> &prop; x</code>, <code>Y<sub>1,-1</sub> &prop; y</code>, <code>Y<sub>1,0</sub> &prop; z</code>.</li>
<li><a href=\"modelica://GeneralRelativity.Harmonics.clebschGordan\">clebschGordan</a> and
<a href=\"modelica://GeneralRelativity.Harmonics.wigner3j\">wigner3j</a> take <code>Real</code>
angular momenta so half-integers can be passed directly (<code>0.5</code>, <code>1.5</code>);
each argument is doubled and rounded to an integer, so pass exact half-integers.
They use the Condon-Shortley phase (<code>&lang;j<sub>1</sub> j<sub>1</sub> j<sub>2</sub> (J - j<sub>1</sub>)|J J&rang; &gt; 0</code>)
and return zero whenever a selection rule fails.</li>
</ul>
<h4>Numerics</h4>
<p>Spherical harmonics use the stable recurrence for the fully normalized functions
<a href=\"modelica://GeneralRelativity.Harmonics.normalizedLegendre\">normalizedLegendre</a>,
so they stay finite for <code>l</code> in the hundreds. The unnormalized
<code>associatedLegendre</code> grows like <code>(2l - 1)!!</code> and overflows a double
near <code>l = 150</code>. Coupling coefficients evaluate the Racah sum with
<a href=\"modelica://GeneralRelativity.Harmonics.logFactorial\">logFactorial</a>; the alternating
sum keeps about twelve significant digits for <code>j &le; 20</code>.
<a href=\"modelica://GeneralRelativity.Harmonics.gaussLegendre\">gaussLegendre</a> nodes in
<code>cos&theta;</code> with uniform <code>&phi;</code> points integrate band-limited functions on
the sphere exactly.</p>
<h4>Use in gravity</h4>
<p>Multipole expansions of the Newtonian potential
(<a href=\"modelica://GeneralRelativity.Perturbation.Potentials.oblateJ2\">J2</a> is the
<code>l = 2, m = 0</code> term), tensor spherical harmonics of black-hole perturbation
theory and the angular integrals of second-order perturbations (Gaunt coefficients)
are built from these functions.</p>
</html>"));
end Harmonics;
