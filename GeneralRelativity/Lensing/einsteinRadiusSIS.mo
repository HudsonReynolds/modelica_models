within GeneralRelativity.Lensing;
function einsteinRadiusSIS "Einstein radius 4 pi (sigma/c)^2 Dls/Ds of a singular isothermal sphere"
  input Real sigma(unit = "m/s") "Line-of-sight velocity dispersion of the lens galaxy";
  input Real Ds(unit = "m") "Distance to the source";
  input Real Dls(unit = "m") "Distance from the lens to the source";
  output Real thetaE(unit = "rad");
algorithm
  thetaE := 4.0 * GeneralRelativity.Constants.pi * (sigma / GeneralRelativity.Constants.c)^2 * Dls / Ds;
  annotation(Inline = true, Documentation(info="<html>
<p>The isothermal sphere <code>&rho; = &sigma;&sup2;/(2&pi;G r&sup2;)</code> has projected mass
<code>M(&lt;R) = &pi;&sigma;&sup2;R/G</code> inside radius <code>R</code>, a flat rotation curve and a constant
deflection <code>&theta;<sub>E</sub></code> for every image position: a good model of a massive galaxy.</p>
</html>"));
end einsteinRadiusSIS;
