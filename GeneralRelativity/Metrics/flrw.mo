within GeneralRelativity.Metrics;
function flrw "Spatially flat FLRW metric in comoving Cartesian coordinates {ct, x, y, z}; p = {ct0, n} or {ct0, n, H}"
  extends GeneralRelativity.Interfaces.MetricFunction;
protected
  Real ct0 = p[1] "Reference time c*t0 at which a = 1";
  Real n = p[2] "Power-law exponent: 2/3 dust, 1/2 radiation";
  Real H "Hubble rate H/c in 1/length; nonzero selects de Sitter expansion";
  Real a "Scale factor";
algorithm
  H := 0.0;
  if size(p, 1) >= 3 then
    H := p[3];
  end if;
  a := if abs(H) > 0.0 then exp(H * (x[1] - ct0)) else (x[1] / ct0)^n;
  g := diagonal({-1.0, a^2, a^2, a^2});
  annotation(Documentation(info="<html>
<p><code>ds&sup2; = -c&sup2;dt&sup2; + a(t)&sup2;(dx&sup2; + dy&sup2; + dz&sup2;)</code> with
cosmic time <code>x[1] = ct</code>. Two scale factors are available:</p>
<ul>
<li><code>p = {ct0, n}</code>: power law <code>a = (ct/ct0)<sup>n</sup></code>, defined for
<code>ct &gt; 0</code>. The Hubble rate is <code>H/c = n/(ct)</code>.
<code>n = 2/3</code> is the Einstein-de Sitter dust universe with
<code>&rho; = 3H&sup2;/(8&pi;G)</code>, <code>n = 1/2</code> is radiation (traceless
Ricci, <code>R = 0</code>). In general the source is a perfect fluid with
<code>w = p/(&rho;c&sup2;) = 2/(3n) - 1</code>.</li>
<li><code>p = {ct0, n, H}</code> with <code>H &ne; 0</code>: de Sitter expansion
<code>a = exp(H(ct - ct0))</code>, <code>H</code> in 1/length; then
<code>R<sub>ab</sub> = &Lambda;g<sub>ab</sub></code> with <code>&Lambda; = 3H&sup2;</code>
and <code>n</code> is ignored.</li>
</ul>
<p>The comoving observer has <code>u = (c, 0, 0, 0)</code>, and the Friedmann
equations read <code>G<sub>00</sub> = 3H&sup2;/c&sup2;</code>,
<code>G<sub>ij</sub> = -(2&auml;/a + H&sup2;) a&sup2;&delta;<sub>ij</sub>/c&sup2;</code>.
For curved spatial sections see
<a href=\"modelica://GeneralRelativity.Metrics.flrwCurved\">flrwCurved</a>.</p>
</html>"));
end flrw;
