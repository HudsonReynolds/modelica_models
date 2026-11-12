within GeneralRelativity.Perturbation;
package Potentials "Analytic Newtonian potentials, gradients and Hessians of weak-field sources"
  annotation(Documentation(info="<html>
<p>Potentials are negative near mass and the acceleration is minus the gradient.
Each potential has a matching gradient with the same parameter vector:</p>
<table>
<tr><th>Source</th><th><code>p</code></th><th>Functions</th></tr>
<tr><td>Point mass at the origin</td><td>{GM}</td>
<td><a href=\"modelica://GeneralRelativity.Perturbation.Potentials.pointMass\">pointMass</a>,
<a href=\"modelica://GeneralRelativity.Perturbation.Potentials.pointMassGradient\">pointMassGradient</a>,
<a href=\"modelica://GeneralRelativity.Perturbation.Potentials.pointMassHessian\">pointMassHessian</a></td></tr>
<tr><td>Oblate body, J2 about +z</td><td>{GM, R, J2}</td>
<td><a href=\"modelica://GeneralRelativity.Perturbation.Potentials.oblateJ2\">oblateJ2</a>,
<a href=\"modelica://GeneralRelativity.Perturbation.Potentials.oblateJ2Gradient\">oblateJ2Gradient</a></td></tr>
<tr><td>Two point masses, second at <code>x2</code></td><td>{GM1, GM2, x2, y2, z2}</td>
<td><a href=\"modelica://GeneralRelativity.Perturbation.Potentials.twoBody\">twoBody</a>,
<a href=\"modelica://GeneralRelativity.Perturbation.Potentials.twoBodyGradient\">twoBodyGradient</a></td></tr>
</table>
<p><a href=\"modelica://GeneralRelativity.Perturbation.Potentials.hessianFromGradient\">hessianFromGradient</a>
differences any gradient function once, which keeps about ten significant
figures and feeds
<a href=\"modelica://GeneralRelativity.Perturbation.weakFieldCurvature\">weakFieldCurvature</a>.</p>
</html>"));
end Potentials;
