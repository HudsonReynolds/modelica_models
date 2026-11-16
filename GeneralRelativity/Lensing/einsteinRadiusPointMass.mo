within GeneralRelativity.Lensing;
function einsteinRadiusPointMass "Einstein radius sqrt(4 G M / c^2 Dls / (Dl Ds)) of a point mass"
  input Real GM(unit = "m3/s2") "Gravitational parameter of the lens";
  input Real Dl(unit = "m") "Distance to the lens";
  input Real Ds(unit = "m") "Distance to the source";
  input Real Dls(unit = "m") "Distance from the lens to the source";
  output Real thetaE(unit = "rad");
algorithm
  thetaE := sqrt(4.0 * GM / GeneralRelativity.Constants.c^2 * Dls / (Dl * Ds));
  annotation(Inline = true, Documentation(info="<html>
<p>A source directly behind the lens is imaged into a ring of this angular radius. For
<code>10<sup>12</sup> M<sub>&#9737;</sub></code> at <code>z<sub>l</sub> = 0.5</code>, <code>z<sub>s</sub> = 2</code>
(flat, <code>H<sub>0</sub> = 70</code>, <code>&Omega;<sub>m</sub> = 0.3</code>) it is 2.03 arcsec.</p>
</html>"));
end einsteinRadiusPointMass;
