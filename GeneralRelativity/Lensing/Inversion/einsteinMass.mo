within GeneralRelativity.Lensing.Inversion;
function einsteinMass "Einstein radius and projected mass inside the critical curve of a lens parameter vector, with standard errors from its covariance"
  input Real lens[7] "Lens parameters (Inversion package documentation)";
  input Real lensCovariance[7, 7] = zeros(7, 7) "Covariance of the lens parameters (fitLens); zeros give zero errors";
  input Boolean pointMass = false;
  input Real Dl(unit = "m") "Distance to the lens";
  input Real Ds(unit = "m") "Distance to the source";
  input Real Dls(unit = "m") "Distance from the lens to the source";
  output Real thetaE(unit = "rad") "Equal-area Einstein radius: b sqrt(q) for the isothermal ellipse, b for a point mass";
  output Real thetaESigma(unit = "rad") "Its standard error";
  output Real GM(unit = "m3/s2") "G times the projected mass inside the critical curve";
  output Real GMSigma(unit = "m3/s2") "Its standard error";
protected
  Real b = max(abs(lens[1]), 1.0e-300);
  Real q = if pointMass then 1.0 else min(max(lens[4], 1.0e-3), 1.0);
  Real relativeVariance "Variance of thetaE / thetaE";
algorithm
  thetaE := b * sqrt(q);
  GM := GeneralRelativity.Lensing.massInEinsteinRadius(thetaE, Dl, Ds, Dls);
  // d(thetaE)/thetaE = db/b + dq/(2q), so its variance is the quadratic form of (1/b, 1/(2q)) with the (b, q) block.
  relativeVariance := lensCovariance[1, 1] / b^2
    + (if pointMass then 0.0 else lensCovariance[4, 4] / (4.0 * q^2) + lensCovariance[1, 4] / (b * q));
  thetaESigma := thetaE * sqrt(max(relativeVariance, 0.0));
  GMSigma := 2.0 * GM * sqrt(max(relativeVariance, 0.0));
  annotation(Documentation(info="<html>
<p>The isothermal ellipse of <a href=\"modelica://GeneralRelativity.Lensing.deflectionSIE\">deflectionSIE</a> has convergence
<code>&kappa; = q b/(2&psi;)</code>, <code>&psi; = &radic;(q&sup2;x&sup2; + y&sup2;)</code>. Inside the iso-density ellipse
<code>&psi; = &psi;<sub>0</sub></code> (area <code>&pi;&psi;<sub>0</sub>&sup2;/q</code>) the projected mass is <code>&pi;b&psi;<sub>0</sub></code>
in units of <code>&Sigma;<sub>cr</sub>D<sub>l</sub>&sup2;</code>, so the mean convergence is 1 inside <code>&psi;<sub>0</sub> = qb</code>: the
critical curve, with semi-axes <code>b</code> and <code>qb</code>, equal-area radius <code>&theta;<sub>E</sub> = b&radic;q</code> and mass
<code>&pi;qb&sup2;&Sigma;<sub>cr</sub>D<sub>l</sub>&sup2;</code>, which is
<a href=\"modelica://GeneralRelativity.Lensing.massInEinsteinRadius\">massInEinsteinRadius</a>(&theta;<sub>E</sub>). That equal-area radius is
also what a circle fitted to the arcs measures
(<a href=\"modelica://GeneralRelativity.Lensing.Inversion.einsteinRadiusFromImage\">einsteinRadiusFromImage</a>), so the radius-based
estimate and the fit are directly comparable. External shear adds no mass. The mass inside a larger aperture <code>R</code> follows
from the isothermal profile: <code>M(R) = GM/G &middot; R/(D<sub>l</sub>&theta;<sub>E</sub>)</code>.</p>
<p>The errors use the covariance of <code>b</code> and <code>q</code> together (they are anti-correlated: a rounder lens needs a smaller
<code>b</code> for the same ring), which <a href=\"modelica://GeneralRelativity.Lensing.Inversion.fitLens\">fitLens</a> returns.</p>
</html>"));
end einsteinMass;
