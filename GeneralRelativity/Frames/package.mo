within GeneralRelativity;
package Frames "Orthonormal and null tetrads, Newman-Penrose Weyl scalars, Petrov classification and energy conditions"
  type PetrovType = enumeration(
    I "Type I: four distinct principal null directions",
    II "Type II: one double principal null direction",
    D "Type D: two double principal null directions",
    III "Type III: one triple principal null direction",
    N "Type N: one quadruple principal null direction",
    O "Type O: Weyl tensor vanishes (conformally flat)") "Petrov type of the Weyl tensor at an event";
  annotation(Documentation(info="<html>
<p>A <em>tetrad</em> is stored as a 4&times;4 array <code>e</code> whose row <code>e[a, :]</code> holds the
contravariant coordinate components <code>e<sub>a</sub><sup>&mu;</sup></code> of frame vector <code>a</code>
(<code>a = 1</code> is the timelike leg <code>e<sub>0</sub></code>). Orthonormality is
<code>e g e<sup>T</sup> = &eta; = diag(&minus;1, 1, 1, 1)</code>, and frame components of a covariant tensor are
<code>T<sub>ab</sub> = e<sub>a</sub><sup>&mu;</sup>e<sub>b</sub><sup>&nu;</sup>T<sub>&mu;&nu;</sub></code>
(<a href=\"modelica://GeneralRelativity.Frames.frameComponents2\">frameComponents2</a>,
<a href=\"modelica://GeneralRelativity.Frames.frameComponents4\">frameComponents4</a>).</p>
<ul>
<li><a href=\"modelica://GeneralRelativity.Frames.orthonormalTetrad\">orthonormalTetrad</a>: the Eulerian
(zero-angular-momentum) observer <code>n<sup>&mu;</sup> = &minus;g<sup>&mu;0</sup>/&radic;(&minus;g<sup>00</sup>)</code>
followed by Gram-Schmidt of <code>&part;<sub>1</sub>, &part;<sub>2</sub>, &part;<sub>3</sub></code>; valid wherever
<code>t = const</code> is spacelike (<code>g<sup>00</sup> &lt; 0</code>), including inside an ergosphere.</li>
<li><a href=\"modelica://GeneralRelativity.Frames.observerTetrad\">observerTetrad</a>: the same with a given
four-velocity as <code>e<sub>0</sub></code>.</li>
<li><a href=\"modelica://GeneralRelativity.Frames.nullTetrad\">nullTetrad</a>:
<code>l = (e<sub>0</sub>+e<sub>1</sub>)/&radic;2</code>, <code>n = (e<sub>0</sub>&minus;e<sub>1</sub>)/&radic;2</code>,
<code>m = (e<sub>2</sub> + i e<sub>3</sub>)/&radic;2</code>, so <code>l&middot;n = &minus;1</code> and
<code>m&middot;m&#772; = 1</code>.</li>
<li><a href=\"modelica://GeneralRelativity.Frames.weylScalars\">weylScalars</a> and
<a href=\"modelica://GeneralRelativity.Frames.weylInvariants\">weylInvariants</a>: Newman-Penrose
<code>&Psi;<sub>0</sub>..&Psi;<sub>4</sub></code> and the invariants <code>I</code>, <code>J</code>.</li>
<li><a href=\"modelica://GeneralRelativity.Frames.weylElectricMagnetic\">weylElectricMagnetic</a>,
<a href=\"modelica://GeneralRelativity.Frames.classifyPetrovMatrix\">classifyPetrovMatrix</a>,
<a href=\"modelica://GeneralRelativity.Frames.classifyPetrov\">classifyPetrov</a> and
<a href=\"modelica://GeneralRelativity.Frames.classifyPetrovCurvature\">classifyPetrovCurvature</a>:
Petrov type (<a href=\"modelica://GeneralRelativity.Frames.PetrovType\">PetrovType</a>) from the eigenstructure of
<code>Q = E + iB</code>.</li>
<li><a href=\"modelica://GeneralRelativity.Frames.energyConditions\">energyConditions</a>: rest-frame energy
density, principal pressures and the null, weak, strong and dominant energy conditions of a
Hawking-Ellis type I stress-energy tensor.</li>
</ul>
<p>Complex quantities are returned as separate real and imaginary arrays so that the package needs
nothing beyond plain Modelica.</p>
</html>"));
end Frames;
