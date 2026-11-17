within GeneralRelativity.CompactBinaries;
function ringdownQuality "Quality factor of the fundamental l = m = 2 quasi-normal mode of a Kerr black hole"
  input Real chi "Dimensionless spin a/M, 0 <= chi < 1";
  output Real Q "omega tau / 2";
algorithm
  Q := 0.7 + 1.4187 * (1.0 - chi)^(-0.4990);
  annotation(Documentation(info="<html>
<p>Fit of Berti, Cardoso and Will (2006): <code>Q = 0.7 + 1.4187(1 &minus; &chi;)<sup>&minus;0.4990</sup></code>. The
amplitude decays as <code>exp(&minus;t/&tau;)</code> with <code>&tau; = 2Q/&omega;</code>: about 11.5 M for
Schwarzschild (exact 11.24 M), 12.2 M<sub>f</sub> at the equal-mass remnant spin 0.686.</p>
</html>"));
end ringdownQuality;
