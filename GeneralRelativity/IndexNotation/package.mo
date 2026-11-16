within GeneralRelativity;
package IndexNotation "Inline tensor expressions in index notation, such as \"g_{mu nu} g^{mu nu}\", evaluated against a system of named tensors at one event"
  annotation(Documentation(info="<html>
<p>This package evaluates tensor expressions written as text in the usual
abstract-index (LaTeX-like) notation against a
<a href=\"modelica://GeneralRelativity.IndexNotation.TensorSystem\">TensorSystem</a>:
a metric plus a table of named tensors at one event. Repeated indices are
summed, indices are raised and lowered automatically with the system's metric,
and the result keeps whatever free indices remain:</p>
<pre>sys = IndexNotation.systemAt(GeneralRelativity.Metrics.schwarzschild, {0, 10, pi/2, 0}, {2.0}, fill(1e-3, 4));
IndexNotation.scalar(sys, \"g_{mu nu} g^{mu nu}\")              // 4
IndexNotation.scalar(sys, \"R_{abcd} R^{abcd}\")                // Kretschmann scalar
IndexNotation.matrix(sys, \"R^{a}_{bac}\")                      // Ricci tensor R_bc
IndexNotation.matrix(sys, \"R_{ab} - 1/2 g_{ab} R - G_{ab}\")   // zero
IndexNotation.scalar(sys, \"epsilon_{abcd} epsilon^{abcd}\")    // -24</pre>

<h4>Expression syntax</h4>
<ul>
<li>A tensor is a name followed by index groups: <code>_</code> lowers,
<code>^</code> raises. Braces are optional: <code>g_{mu nu}</code>,
<code>g_{\\mu\\nu}</code>, <code>g_{ab}</code>, <code>g_ab</code> and
<code>R^a_{bcd}</code> all parse. Inside a group, tokens are separated by
whitespace, commas or backslashes; a token that is a Greek letter name
(<code>mu</code>, <code>nu</code>, <code>alpha</code>, ...) is one index and any
other token is one index per character, so <code>abcd</code> is four indices.</li>
<li>An index repeated twice, in one tensor or across a product, is summed
(contracted). If both copies sit in the same position the metric is inserted,
so <code>R_{ab} R_{ab}</code> equals <code>R_{ab} R^{ab}</code>. An index may
not appear more than twice.</li>
<li>A digit <code>0..3</code> is a fixed coordinate index in physics labels:
<code>g_{00}</code> is the <code>g[1, 1]</code> component and
<code>R_{0101}</code> is <code>riemann[1, 2, 1, 2]</code>; the position
(<code>^</code> or <code>_</code>) is honoured before the component is taken.</li>
<li>Juxtaposition or <code>*</code> multiplies, <code>/</code> divides by a
scalar-valued factor, <code>+</code> and <code>-</code> add terms whose free
indices agree as sets (the order and the positions may differ; the second
term is permuted and raised or lowered to match the first), parentheses group,
and real literals such as <code>2</code>, <code>0.5</code> or <code>1e-3</code>
scale. Multiplication binds tighter than addition.</li>
<li>Names are case sensitive identifiers of letters and digits. The same name may
hold tensors of different ranks, chosen by the number of indices:
<code>R</code> is the Riemann tensor with four indices, the Ricci tensor with
two and the Ricci scalar with none.</li>
</ul>

<h4>Built-in names</h4>
<p><a href=\"modelica://GeneralRelativity.IndexNotation.systemAt\">systemAt</a> and
<a href=\"modelica://GeneralRelativity.IndexNotation.systemFromCurvature\">systemFromCurvature</a>
define, all in the conventions of
<a href=\"modelica://GeneralRelativity.Tensors\">Tensors</a>:</p>
<table>
<tr><th>Name</th><th>Indices</th><th>Meaning</th></tr>
<tr><td><code>g</code>, <code>delta</code></td><td>2</td><td>metric; <code>g^{ab}</code> is the inverse and <code>delta^a_b</code> the identity</td></tr>
<tr><td><code>Gamma</code></td><td>3</td><td>Christoffel symbols <code>Gamma^i_{jk}</code> (stored with the first index lowered, so raise it)</td></tr>
<tr><td><code>R</code>, <code>Riemann</code></td><td>4</td><td><code>R_abcd</code></td></tr>
<tr><td><code>R</code>, <code>Ric</code>, <code>Ricci</code></td><td>2</td><td><code>R_ab</code></td></tr>
<tr><td><code>R</code></td><td>0</td><td>Ricci scalar</td></tr>
<tr><td><code>G</code>, <code>Einstein</code></td><td>2</td><td><code>G_ab</code></td></tr>
<tr><td><code>C</code>, <code>Weyl</code></td><td>4</td><td><code>C_abcd</code></td></tr>
<tr><td><code>P</code>, <code>Schouten</code></td><td>2</td><td>Schouten tensor</td></tr>
<tr><td><code>S</code></td><td>2</td><td>trace-free Ricci tensor</td></tr>
<tr><td><code>T</code></td><td>2</td><td>stress-energy <code>T_ab</code> implied by the Einstein tensor for the given <code>G</code>, <code>c</code>, <code>Lambda</code></td></tr>
<tr><td><code>K</code></td><td>0</td><td>Kretschmann scalar</td></tr>
<tr><td><code>epsilon</code>, <code>eps</code></td><td>4</td><td>Levi-Civita tensor <code>sqrt(|g|) [abcd]</code></td></tr>
</table>
<p><a href=\"modelica://GeneralRelativity.IndexNotation.systemFromMetric\">systemFromMetric</a>
defines only <code>g</code>, <code>delta</code> and <code>epsilon</code>, for
algebra on user tensors. Add your own with
<a href=\"modelica://GeneralRelativity.IndexNotation.defineScalar\">defineScalar</a>,
<a href=\"modelica://GeneralRelativity.IndexNotation.defineVector\">defineVector</a>,
<a href=\"modelica://GeneralRelativity.IndexNotation.defineMatrix\">defineMatrix</a>,
<a href=\"modelica://GeneralRelativity.IndexNotation.defineTensor3\">defineTensor3</a>,
<a href=\"modelica://GeneralRelativity.IndexNotation.defineTensor4\">defineTensor4</a>
or the flat <a href=\"modelica://GeneralRelativity.IndexNotation.define\">define</a>;
each takes the index positions of the components you pass and stores the
covariant form. Reusing a name and rank replaces the earlier tensor.</p>

<h4>Results</h4>
<p><a href=\"modelica://GeneralRelativity.IndexNotation.evaluate\">evaluate</a>
returns a <a href=\"modelica://GeneralRelativity.IndexNotation.TensorValue\">TensorValue</a>
with the free index labels, their positions and the flattened components;
<a href=\"modelica://GeneralRelativity.IndexNotation.scalar\">scalar</a>,
<a href=\"modelica://GeneralRelativity.IndexNotation.vector\">vector</a>,
<a href=\"modelica://GeneralRelativity.IndexNotation.matrix\">matrix</a>,
<a href=\"modelica://GeneralRelativity.IndexNotation.tensor3\">tensor3</a> and
<a href=\"modelica://GeneralRelativity.IndexNotation.tensor4\">tensor4</a> check
the rank and return ordinary arrays, with the free indices in the order they
first appear in the expression and in the positions the expression gave them.
<a href=\"modelica://GeneralRelativity.IndexNotation.toString\">toString</a>
formats a value for printing. A syntax error, an unknown name, an index used
three times or a mismatch of free indices between added terms stops with an
assertion that quotes the expression and the character position.</p>

<h4>Limits</h4>
<p>The evaluator is numeric, in four dimensions and at one event: it is a
calculator for components, not a symbolic tensor algebra, and it knows no
derivatives (take <code>Gamma</code> and the curvature tensors from the system
instead). Intermediate values may carry at most
<code>Internal.MaxRank = 6</code> free indices, so <code>R_{abcd} R_{efgh}</code>
without a contraction is refused. Names and components are held in fixed-size
tables sized by the constants in
<a href=\"modelica://GeneralRelativity.IndexNotation.Internal\">Internal</a>.</p>
</html>"));
end IndexNotation;
