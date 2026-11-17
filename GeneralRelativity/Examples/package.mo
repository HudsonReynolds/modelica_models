within GeneralRelativity;
package Examples "Worked models: Schwarzschild orbit and curvature, Newton-Cartan comparison, Earth J2, cislunar, Poisson solvers, star interior, wave pulse"
  annotation(Documentation(info="<html>
<table>
<tr><th>Model</th><th>Shows</th></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Examples.SchwarzschildOrbit\">SchwarzschildOrbit</a></td>
<td>Geodesic integration, Killing constants, curvature at the particle.</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Examples.SchwarzschildNewtonCartanOrbit\">SchwarzschildNewtonCartanOrbit</a></td>
<td>The same orbit beside a Newton-Cartan geodesic from identical initial data, advanced in the same coordinate time: precession against a closed Kepler ellipse.</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Examples.SchwarzschildOrbitScene\">SchwarzschildOrbitScene</a>,
<a href=\"modelica://GeneralRelativity.Examples.EarthJ2OrbitScene\">EarthJ2OrbitScene</a>,
<a href=\"modelica://GeneralRelativity.Examples.CislunarOrbitScene\">CislunarOrbitScene</a></td>
<td>The same orbits plus the Three.js scene stream (need <code>Modelica.Utilities</code> to link); every scene can be viewed as embedding, river, clocks, spacetime diagram or (with an analytic metric) ray tracer.</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Examples.KerrOrbitScene\">KerrOrbitScene</a></td>
<td>Equatorial Kerr geodesic with a SceneRecorder: frame-dragged river, ray-traced asymmetric shadow.</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Examples.TwinParadox\">TwinParadox</a></td>
<td>Twin paradox in Minkowski space with constant proper acceleration; ages checked against hyperbolic motion; interactive diagram in the viewer.</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Examples.LightDeflection\">LightDeflection</a></td>
<td>A Schwarzschild null geodesic's bending angle against 4GM/(c&sup2;b) and its higher orders.</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Examples.GalaxyLensing\">GalaxyLensing</a></td>
<td>Galaxy-scale strong lensing: distances, Einstein radii, a microlensing light curve; interactive lensed sky in the viewer.</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Examples.LensInversion\">LensInversion</a>,
<a href=\"modelica://GeneralRelativity.Examples.LensInversionFromFile\">LensInversionFromFile</a></td>
<td>The inverse: from a lensed image (synthetic with a known mass, or read from a file) to the lens mass: ring radius, mass estimate, simulation, &chi;&sup2; and a least-squares fit with standard errors; observed, simulated and residual images in the viewer.</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Examples.SchwarzschildCurvature\">SchwarzschildCurvature</a></td>
<td>Christoffel symbols against closed form, vacuum Ricci, Kretschmann, implied stress-energy.</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.Examples.EarthJ2Orbit\">EarthJ2Orbit</a></td>
<td>Weak-field Earth: J2, post-Newtonian and Lense-Thirring accelerations, tidal tensor.</td></tr>
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
