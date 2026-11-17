within GeneralRelativity.Lensing;
function massInEinsteinRadius "Projected mass inside an Einstein radius, pi (Dl thetaE)^2 Sigma_cr, as a gravitational parameter: the inverse of einsteinRadiusPointMass"
  input Real thetaE(unit = "rad") "Einstein radius";
  input Real Dl(unit = "m") "Distance to the lens";
  input Real Ds(unit = "m") "Distance to the source";
  input Real Dls(unit = "m") "Distance from the lens to the source";
  output Real GM(unit = "m3/s2") "G times the projected mass inside thetaE";
algorithm
  GM := GeneralRelativity.Constants.c^2 * thetaE^2 * Dl * Ds / (4.0 * Dls);
  annotation(Inline = true, Documentation(info="<html>
<p>Any axially symmetric lens images a source on its axis into a ring whose radius encloses a mean
convergence of 1, so the projected mass inside the ring is
<code>M = &pi; (D<sub>l</sub>&theta;<sub>E</sub>)&sup2; &Sigma;<sub>cr</sub> = c&sup2;&theta;<sub>E</sub>&sup2;D<sub>l</sub>D<sub>s</sub>/(4G D<sub>ls</sub>)</code>
whatever the mass profile. For a point mass it is the whole mass
(<a href=\"modelica://GeneralRelativity.Lensing.einsteinRadiusPointMass\">einsteinRadiusPointMass</a> inverted); for
an isothermal sphere <a href=\"modelica://GeneralRelativity.Lensing.velocityDispersionSIS\">velocityDispersionSIS</a>
gives &sigma; and the mass inside any other radius <code>R</code> follows from <code>M(&lt;R) = &pi;&sigma;&sup2;R/G</code>.
This is the first, model-independent mass estimate of
<a href=\"modelica://GeneralRelativity.Lensing.Inversion\">Lensing.Inversion</a>.</p>
</html>"));
end massInEinsteinRadius;
