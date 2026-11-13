within GeneralRelativity;
package ExteriorCalculus "Differential forms in components: wedge, Hodge dual, interior product, d, codifferential and Cartan's structure equations"
  annotation(Documentation(info="<html>
<p>Four-dimensional differential forms with lower coordinate indices, stored as
fully antisymmetric component arrays: a 1-form <code>w[4]</code>, a 2-form
<code>F[4,4]</code>, a 3-form <code>H[4,4,4]</code> and a 4-form <code>V[4,4,4,4]</code>.</p>

<h4>Normalization</h4>
<ul>
<li>Components: <code>&omega; = (1/p!) &omega;<sub>a<sub>1</sub>..a<sub>p</sub></sub>
dx<sup>a<sub>1</sub></sup>&and;..&and;dx<sup>a<sub>p</sub></sup></code>.</li>
<li>Wedge: <code>(&alpha;&and;&beta;)<sub>a<sub>1</sub>..a<sub>p</sub>b<sub>1</sub>..b<sub>q</sub></sub>
= ((p+q)!/(p!q!)) &alpha;<sub>[a<sub>1</sub>..a<sub>p</sub></sub>&beta;<sub>b<sub>1</sub>..b<sub>q</sub>]</sub></code>,
so <code>(a&and;b)<sub>ij</sub> = a<sub>i</sub>b<sub>j</sub> - a<sub>j</sub>b<sub>i</sub></code>
(Wald, Carroll, MTW).</li>
<li>Exterior derivative: <code>(d&omega;)<sub>a<sub>0</sub>..a<sub>p</sub></sub> = (p+1) &part;<sub>[a<sub>0</sub></sub>&omega;<sub>a<sub>1</sub>..a<sub>p</sub>]</sub></code>,
e.g. <code>(dA)<sub>ab</sub> = &part;<sub>a</sub>A<sub>b</sub> - &part;<sub>b</sub>A<sub>a</sub></code>.
Field derivatives are central finite differences.</li>
<li>Interior product: <code>(i<sub>X</sub>&omega;)<sub>a<sub>2</sub>..</sub> = X<sup>a</sup>&omega;<sub>a a<sub>2</sub>..</sub></code>, an antiderivation.</li>
<li>Levi-Civita tensor: <code>&epsilon;<sub>abcd</sub> = &radic;|det g| [abcd]</code> with
<code>[0123] = +1</code> (Modelica indices <code>[1234]</code>), so
<code>&epsilon;<sup>0123</sup> = -1/&radic;|g|</code>.</li>
<li>Hodge dual: <code>(*&omega;)<sub>b<sub>1</sub>..b<sub>4-p</sub></sub> = (1/p!)
&omega;<sup>a<sub>1</sub>..a<sub>p</sub></sup> &epsilon;<sub>a<sub>1</sub>..a<sub>p</sub>b<sub>1</sub>..b<sub>4-p</sub></sub></code>
(Carroll 2.9). Then <code>&alpha;&and;*&beta; = &lang;&alpha;,&beta;&rang; &epsilon;</code> with
<code>&lang;&alpha;,&beta;&rang; = (1/p!)&alpha;<sub>a..</sub>&beta;<sup>a..</sup></code>, and in
Lorentzian signature <code>** = (-1)<sup>p(4-p)+1</sup></code>: <code>-1, +1, -1, +1, -1</code>
for <code>p = 0..4</code>. <code>*1 = &epsilon;</code>, <code>*&epsilon; = -1</code>.</li>
<li>Codifferential: <code>&delta; = -div</code>, i.e. <code>&delta;w = -&nabla;<sup>a</sup>w<sub>a</sub></code>,
<code>(&delta;F)<sub>b</sub> = -&nabla;<sup>a</sup>F<sub>ab</sub></code>; with these
conventions <code>&delta; = *d*</code> in four Lorentzian dimensions.</li>
</ul>

<h4>Maxwell</h4>
<p>With <code>F = dA</code> the Maxwell equations are <code>dF = 0</code> and
<code>&delta;F = &plusmn;&mu;<sub>0</sub>J</code>; the Coulomb field of a charge on
the Reissner-Nordstr&ouml;m metric, <code>F = (Q/r&sup2;) dx<sup>0</sup>&and;dr</code>,
satisfies both with <code>J = 0</code>.</p>

<h4>Tetrads and Cartan's structure equations</h4>
<p>A <a href=\"modelica://GeneralRelativity.ExteriorCalculus.TetradField\">TetradField</a>
returns <code>e[a, &mu;] = e<sub>a</sub><sup>&mu;</sup></code> with
<code>g(e<sub>a</sub>, e<sub>b</sub>) = &eta;<sub>ab</sub> = diag(-1, 1, 1, 1)</code>.
The coframe is <code>&theta;<sup>a</sup> = &eta;<sup>ab</sup>g(e<sub>b</sub>, &middot;)</code>,
the spin connection is defined by <code>&nabla;e<sub>b</sub> = e<sub>a</sub> &otimes; &omega;<sup>a</sup><sub>b</sub></code>,
and the structure equations are</p>
<pre>d&theta;<sup>a</sup> + &omega;<sup>a</sup><sub>b</sub>&and;&theta;<sup>b</sup> = 0,
&Omega;<sup>a</sup><sub>b</sub> = d&omega;<sup>a</sup><sub>b</sub> + &omega;<sup>a</sup><sub>c</sub>&and;&omega;<sup>c</sup><sub>b</sub> = &frac12;R<sup>a</sup><sub>b&mu;&nu;</sub>dx<sup>&mu;</sup>&and;dx<sup>&nu;</sup>,</pre>
<p>with <code>R<sup>a</sup><sub>b&mu;&nu;</sub> = &theta;<sup>a</sup><sub>&rho;</sub>R<sup>&rho;</sup><sub>&sigma;&mu;&nu;</sub>e<sub>b</sub><sup>&sigma;</sup></code>
in the Riemann sign convention of the library.
<a href=\"modelica://GeneralRelativity.ExteriorCalculus.cartanTorsion\">cartanTorsion</a> and
<a href=\"modelica://GeneralRelativity.ExteriorCalculus.cartanCurvatureResidual\">cartanCurvatureResidual</a>
return the residuals of both equations;
<a href=\"modelica://GeneralRelativity.ExteriorCalculus.schwarzschildTetrad\">schwarzschildTetrad</a>
is the static frame of Schwarzschild.</p>

<h4>Field signatures</h4>
<p><a href=\"modelica://GeneralRelativity.ExteriorCalculus.ScalarField\">ScalarField</a>,
<a href=\"modelica://GeneralRelativity.ExteriorCalculus.OneFormField\">OneFormField</a>,
<a href=\"modelica://GeneralRelativity.ExteriorCalculus.TwoFormField\">TwoFormField</a>,
<a href=\"modelica://GeneralRelativity.ExteriorCalculus.ThreeFormField\">ThreeFormField</a> and
<a href=\"modelica://GeneralRelativity.ExteriorCalculus.TetradField\">TetradField</a> take
<code>(x[4], p[:])</code> like a metric and are passed by name to the
differential operators.</p>
</html>"));
end ExteriorCalculus;
