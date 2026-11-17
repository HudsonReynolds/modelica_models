within GeneralRelativity.CompactBinaries;
function ringdownOmega "Angular frequency of the fundamental l = m = 2 quasi-normal mode of a Kerr black hole, in units of 1/M"
  input Real chi "Dimensionless spin a/M, 0 <= chi < 1";
  output Real omega "M omega";
algorithm
  omega := 1.5251 - 1.1568 * (1.0 - chi)^0.1292;
  annotation(Documentation(info="<html>
<p>Fit of E. Berti, V. Cardoso and C. M. Will, <i>Phys. Rev. D</i> 73, 064030 (2006), table VIII:
<code>M&omega; = 1.5251 &minus; 1.1568(1 &minus; &chi;)<sup>0.1292</sup></code> (within a few percent; 0.368 for
Schwarzschild against the exact 0.3737). The ringing wave has angular frequency <code>&omega;</code>; the
damping time is <code>2Q/&omega;</code>, see <a href=\"modelica://GeneralRelativity.CompactBinaries.ringdownQuality\">ringdownQuality</a>.</p>
</html>"));
end ringdownOmega;
