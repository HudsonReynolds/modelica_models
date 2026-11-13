within GeneralRelativity;
package Foliation "3+1 (ADM) split: lapse, shift, spatial metric, extrinsic curvature, 3-curvature and the Hamiltonian and momentum constraints"
  annotation(Documentation(info="<html>
<p>A spacetime chart whose time coordinate <code>x[1] = ct</code> has spacelike
level sets is foliated into slices <code>&Sigma;<sub>t</sub></code>. The metric
is written in Arnowitt-Deser-Misner form</p>
<pre>ds&sup2; = -&alpha;&sup2; dx<sup>0</sup>dx<sup>0</sup> + &gamma;<sub>ij</sub>(dx<sup>i</sup> + &beta;<sup>i</sup>dx<sup>0</sup>)(dx<sup>j</sup> + &beta;<sup>j</sup>dx<sup>0</sup>),</pre>
<p>so <code>g<sub>00</sub> = -&alpha;&sup2; + &beta;<sub>k</sub>&beta;<sup>k</sup></code>,
<code>g<sub>0i</sub> = &beta;<sub>i</sub></code> and <code>g<sub>ij</sub> = &gamma;<sub>ij</sub></code>,
with lapse <code>&alpha;</code>, shift <code>&beta;<sup>i</sup></code>
(lowered with <code>&gamma;</code>) and spatial metric <code>&gamma;<sub>ij</sub></code>.
Spatial indices <code>i, j = 1..3</code> are the Modelica array indices of 3-vectors and
correspond to spacetime indices <code>2..4</code>. The future unit normal is
<code>n<sub>a</sub> = (-&alpha;, 0, 0, 0)</code>,
<code>n<sup>a</sup> = (1, -&beta;<sup>i</sup>)/&alpha;</code>.</p>

<h4>Sign of the extrinsic curvature</h4>
<p>This library follows MTW, York and Alcubierre:</p>
<pre>K<sub>ij</sub> = -&gamma;<sub>i</sub><sup>a</sup>&gamma;<sub>j</sub><sup>b</sup>&nabla;<sub>a</sub>n<sub>b</sub>
       = -(1/(2&alpha;)) (&part;<sub>0</sub>&gamma;<sub>ij</sub> - D<sub>i</sub>&beta;<sub>j</sub> - D<sub>j</sub>&beta;<sub>i</sub>),   &part;<sub>0</sub> = &part;/&part;(ct),</pre>
<p>so <code>K = &gamma;<sup>ij</sup>K<sub>ij</sub> = -&nabla;<sub>a</sub>n<sup>a</sup></code>:
slices of an expanding universe have <code>K &lt; 0</code> (flat FLRW gives
<code>K = -3H/c</code>) and slices whose normal observers converge have
<code>K &gt; 0</code>. In coordinates, <code>K<sub>ij</sub> = -&alpha;&Gamma;<sup>0</sup><sub>ij</sub>
= n<sup>l</sup>&Gamma;<sub>lij</sub></code>.</p>

<h4>Gauss-Codazzi constraints</h4>
<p>With <code>D</code> the Levi-Civita derivative of <code>&gamma;</code>,
<code>R<sup>(3)</sup></code> its Ricci scalar and <code>G<sub>ab</sub></code> the
4-dimensional Einstein tensor, the contracted Gauss and Codazzi equations are</p>
<pre>H   = R<sup>(3)</sup> + K&sup2; - K<sub>ij</sub>K<sup>ij</sup> - 2 G<sub>ab</sub>n<sup>a</sup>n<sup>b</sup>       = 0,
M<sub>i</sub> = D<sub>j</sub>K<sup>j</sup><sub>i</sub> - D<sub>i</sub>K + G<sub>ab</sub>n<sup>a</sup>&delta;<sup>b</sup><sub>i</sub>               = 0.</pre>
<p>They hold for <i>every</i> metric; Einstein's equations turn them into the
Hamiltonian and momentum constraints
<code>R<sup>(3)</sup> + K&sup2; - K<sub>ij</sub>K<sup>ij</sup> = 16&pi;G&rho;/c<sup>4</sup></code>,
<code>D<sub>j</sub>K<sup>j</sup><sub>i</sub> - D<sub>i</sub>K = 8&pi;G j<sub>i</sub>/c<sup>4</sup></code>
with Eulerian energy density <code>&rho; = T<sub>ab</sub>n<sup>a</sup>n<sup>b</sup></code>
and momentum density <code>j<sub>i</sub> = -T<sub>ab</sub>n<sup>a</sup>&delta;<sup>b</sup><sub>i</sub></code>.
<a href=\"modelica://GeneralRelativity.Foliation.admAt\">admAt</a> evaluates both
residuals against the 4-dimensional Einstein tensor of the same finite-difference
jet, so they vanish to rounding error for any metric and are a check of the split
itself; the source terms are returned separately.</p>

<h4>Functions</h4>
<ul>
<li><a href=\"modelica://GeneralRelativity.Foliation.admDecomposition\">admDecomposition</a>: lapse, shift and spatial metric (algebra);</li>
<li><a href=\"modelica://GeneralRelativity.Foliation.normalVector\">normalVector</a>: unit normal <code>n<sup>a</sup></code>, <code>n<sub>a</sub></code>;</li>
<li><a href=\"modelica://GeneralRelativity.Foliation.extrinsicCurvature\">extrinsicCurvature</a>: <code>K<sub>ij</sub></code> from <code>g</code> and <code>&part;g</code>;</li>
<li><a href=\"modelica://GeneralRelativity.Foliation.extrinsicCurvatureGradient\">extrinsicCurvatureGradient</a>: <code>&part;<sub>k</sub>K<sub>ij</sub></code> from the 2-jet;</li>
<li><a href=\"modelica://GeneralRelativity.Foliation.admFromJet\">admFromJet</a> and
<a href=\"modelica://GeneralRelativity.Foliation.admAt\">admAt</a>: the full
<a href=\"modelica://GeneralRelativity.Foliation.ADM\">ADM</a> record.</li>
</ul>
<p>References: Arnowitt, Deser and Misner (1962); Misner, Thorne and Wheeler, ch. 21;
Alcubierre, <i>Introduction to 3+1 Numerical Relativity</i>, ch. 2.</p>
</html>"));
end Foliation;
