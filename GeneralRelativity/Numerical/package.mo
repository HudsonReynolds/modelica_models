within GeneralRelativity;
package Numerical "Finite-difference and finite-element field solvers for weak-field and linearized problems"
  annotation(Documentation(info="<html>
<p>The Newtonian limit of the field equations is Poisson's equation
<code>&nabla;&sup2;&Phi; = 4&pi;G&rho;</code>; the linearized equations in the Lorenz
gauge are wave equations. Both are solved here on fixed grids:</p>
<ul>
<li><a href=\"modelica://GeneralRelativity.Numerical.poissonFiniteDifference3D\">poissonFiniteDifference3D</a>:
seven-point stencil with successive over-relaxation on a Cartesian box and Dirichlet
boundary values (for example the monopole far field from
<a href=\"modelica://GeneralRelativity.Numerical.monopoleBoundary\">monopoleBoundary</a>).
Returns the sweep count and residual so the caller can judge convergence.</li>
<li><a href=\"modelica://GeneralRelativity.Numerical.poissonFiniteElementSpherical\">poissonFiniteElementSpherical</a>:
linear finite elements for the spherically symmetric weak form
<code>&int; r&sup2; &Phi;' v' dr = -4&pi;G &int; r&sup2; &rho; v dr</code> with the natural
zero-flux condition at the centre and the exact monopole value at the outer
node.</li>
<li><a href=\"modelica://GeneralRelativity.Numerical.LinearizedWave1D\">LinearizedWave1D</a>:
method-of-lines wave equation for a trace-reversed perturbation component with
Sommerfeld outflow boundaries, integrated by the Modelica solver.</li>
</ul>
<p>The solvers are dimension generic in the array sizes. Loop bounds such as
<code>maxIterations</code> should be constants or evaluated parameters for
compilers that need translation-time bounds.</p>
</html>"));
end Numerical;
