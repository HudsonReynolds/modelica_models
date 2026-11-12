within GeneralRelativity;
package Examples "Worked models: Schwarzschild orbit and curvature, Earth J2, cislunar, Poisson solvers, star interior, wave pulse"
  annotation(Documentation(info="<html>
<table>
<tr><th>Model</th><th>Shows</th></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Examples.SchwarzschildOrbit\">SchwarzschildOrbit</a></td>
<td>Geodesic integration, Killing constants, curvature at the particle, Three.js scene stream.</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Examples.SchwarzschildCurvature\">SchwarzschildCurvature</a></td>
<td>Christoffel symbols against closed form, vacuum Ricci, Kretschmann, implied stress-energy.</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Examples.EarthJ2Orbit\">EarthJ2Orbit</a></td>
<td>Weak-field Earth: J2, post-Newtonian and Lense-Thirring accelerations, tidal field slice.</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Examples.CislunarOrbit\">CislunarOrbit</a></td>
<td>Earth-Moon weak field with a moving Moon, Jacobi constant as the integration check.</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Examples.WeakFieldPoisson\">WeakFieldPoisson</a></td>
<td>Finite-difference and finite-element Poisson solutions of a homogeneous sphere.</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Examples.StarInterior\">StarInterior</a></td>
<td>Tolman-Oppenheimer-Volkoff integration against the interior Schwarzschild solution.</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Examples.GravitationalWavePulse\">GravitationalWavePulse</a></td>
<td>Method-of-lines propagation of a linearized perturbation with outflow boundaries.</td></tr>
</table>
<p>Every example carries continuous assertions so it doubles as a regression
test; <code>Tests/run-general-relativity.mos</code> simulates them.</p>
</html>"));
end Examples;
