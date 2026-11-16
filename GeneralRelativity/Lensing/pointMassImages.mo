within GeneralRelativity.Lensing;
function pointMassImages "Both images of a point-mass lens and their signed magnifications"
  input Real beta(unit = "rad") "Source distance from the lens axis (>= 0)";
  input Real thetaE(unit = "rad") "Einstein radius";
  output Real thetaPlus(unit = "rad") "Outer image, on the source's side";
  output Real thetaMinus(unit = "rad") "Inner image, on the opposite side (negative)";
  output Real muPlus "Magnification of the outer image (> 1)";
  output Real muMinus "Signed magnification of the inner image (< 0: mirror-reversed)";
algorithm
  thetaPlus := 0.5 * (beta + sqrt(beta^2 + 4.0 * thetaE^2));
  thetaMinus := 0.5 * (beta - sqrt(beta^2 + 4.0 * thetaE^2));
  muPlus := 1.0 / (1.0 - (thetaE / thetaPlus)^4);
  muMinus := 1.0 / (1.0 - (thetaE / thetaMinus)^4);
  annotation(Documentation(info="<html>
<p>Solutions of <code>&beta; = &theta; &minus; &theta;<sub>E</sub>&sup2;/&theta;</code> with
<code>&mu; = [1 &minus; (&theta;<sub>E</sub>/&theta;)<sup>4</sup>]<sup>&minus;1</sup></code>. The total brightening
<code>|&mu;<sub>+</sub>| + |&mu;<sub>&minus;</sub>| = (u&sup2; + 2)/(u&radic;(u&sup2; + 4))</code>, <code>u = &beta;/&theta;<sub>E</sub></code>,
is the microlensing light curve; the signed sum is exactly 1.</p>
</html>"));
end pointMassImages;
