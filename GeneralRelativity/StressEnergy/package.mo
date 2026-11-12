within GeneralRelativity;
package StressEnergy "Covariant stress-energy tensors T_ab of standard matter models"
  annotation(Documentation(info="<html>
<p>Every function returns the covariant components <code>T<sub>ab</sub></code> in
the (&minus;+++) signature with <code>x[1] = ct</code>. The four-velocity
<code>u</code> is contravariant with <code>g(u,u) = -c&sup2;</code>; use
<a href=\"modelica://GeneralRelativity.Geodesics.fourVelocity\">Geodesics.fourVelocity</a>
to build it from a coordinate velocity. Energy density measured by the fluid
is <code>T<sub>ab</sub>u<sup>a</sup>u<sup>b</sup>/c&sup2;</code>
(<a href=\"modelica://GeneralRelativity.StressEnergy.energyDensity\">energyDensity</a>).</p>
<table>
<tr><th>Source</th><th>Function</th><th>T<sub>ab</sub></th></tr>
<tr><td>Vacuum</td><td><a href=\"modelica://GeneralRelativity.StressEnergy.vacuum\">vacuum</a></td><td>0</td></tr>
<tr><td>Non-interacting matter</td><td><a href=\"modelica://GeneralRelativity.StressEnergy.dust\">dust</a></td><td>&rho; u<sub>a</sub>u<sub>b</sub></td></tr>
<tr><td>Perfect fluid</td><td><a href=\"modelica://GeneralRelativity.StressEnergy.perfectFluid\">perfectFluid</a></td><td>(&rho; + p/c&sup2;) u<sub>a</sub>u<sub>b</sub> + p g<sub>ab</sub></td></tr>
<tr><td>Radiation</td><td><a href=\"modelica://GeneralRelativity.StressEnergy.radiation\">radiation</a></td><td>perfect fluid with p = &rho;c&sup2;/3</td></tr>
<tr><td>Cosmological constant</td><td><a href=\"modelica://GeneralRelativity.StressEnergy.cosmologicalConstant\">cosmologicalConstant</a></td><td>-&Lambda;c<sup>4</sup>/(8&pi;G) g<sub>ab</sub></td></tr>
</table>
<p><a href=\"modelica://GeneralRelativity.StressEnergy.einsteinResidual\">einsteinResidual</a>
measures how far a metric and a source are from satisfying
<code>G<sub>ab</sub> + &Lambda;g<sub>ab</sub> = (8&pi;G/c<sup>4</sup>) T<sub>ab</sub></code>.</p>
</html>"));
end StressEnergy;
