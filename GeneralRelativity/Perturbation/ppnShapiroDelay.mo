within GeneralRelativity.Perturbation;
function ppnShapiroDelay "One-way Shapiro time delay (1 + gamma) GM/c^3 ln((r1 + r2 + R)/(r1 + r2 - R))"
  input Real GM(unit = "m3/s2") "Gravitational parameter of the deflecting body";
  input Real r1(unit = "m") "Distance of the emitter from the body";
  input Real r2(unit = "m") "Distance of the receiver from the body";
  input Real R(unit = "m") "Emitter-receiver distance";
  input Real gamma = 1.0 "PPN gamma";
  input Real c(unit = "m/s") = GeneralRelativity.Constants.c "Speed of light";
  output Real delay(unit = "s") "Excess coordinate light-travel time over R/c";
algorithm
  delay := (1.0 + gamma) * GM / c^3 * log((r1 + r2 + R) / (r1 + r2 - R));
  annotation(Documentation(info="<html>
<p>Integral of the PPN refractive index <code>1 + (1 + &gamma;)U</code> along the
straight line between emitter and receiver. Near superior conjunction, with
impact parameter <code>b &lt;&lt; r<sub>1</sub>, r<sub>2</sub></code>, it is
<code>(1 + &gamma;)(GM/c&sup3;) ln(4r<sub>1</sub>r<sub>2</sub>/b&sup2;)</code>; a radar
round trip is twice the one-way value (about 250 &micro;s for Earth-Mars grazing the
Sun).</p>
</html>"));
end ppnShapiroDelay;
