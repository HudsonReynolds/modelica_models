within GeneralRelativity.StressEnergy;
function polytropeEnergyDensity "Total energy density of a polytrope, e = rho c^2 + p/(Gamma - 1)"
  input Real rho "Rest-mass density (>= 0)";
  input Real K "Polytropic constant";
  input Real Gamma "Adiabatic index (> 1)";
  input Real c = 1.0 "Speed of light; 1 in geometric units";
  output Real energyDensity "Rest-mass plus internal energy per unit volume";
algorithm
  energyDensity := max(rho, 0.0) * c^2 + K * max(rho, 0.0)^Gamma / (Gamma - 1.0);
  annotation(Documentation(info="<html>
<p>For an adiabatic (isentropic) polytrope the first law
<code>d(e/&rho;) = -p d(1/&rho;)</code> gives the specific internal energy
<code>p/((&Gamma; - 1)&rho;)</code>, hence <code>e = &rho;c&sup2; + p/(&Gamma; - 1)</code>
and <code>de/d&rho; = (e + p)/&rho;</code>. This <code>e/c&sup2;</code> is the
<code>rho</code> argument of
<a href=\"modelica://GeneralRelativity.StressEnergy.perfectFluid\">perfectFluid</a> and
the source of the TOV mass equation in
<a href=\"modelica://GeneralRelativity.Examples.PolytropeStar\">Examples.PolytropeStar</a>.
Some authors instead write the polytrope in the total energy density,
<code>p = K e<sup>&Gamma;</sup></code>; this library uses the rest-mass form above, which
reduces to the Newtonian Lane-Emden polytrope for <code>p &lt;&lt; &rho;c&sup2;</code>.</p>
</html>"));
end polytropeEnergyDensity;
