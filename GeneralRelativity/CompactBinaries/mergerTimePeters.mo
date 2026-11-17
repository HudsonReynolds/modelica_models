within GeneralRelativity.CompactBinaries;
function mergerTimePeters "Time for a circular binary to shrink to zero separation at leading order (Peters 1964)"
  input Real r0 "Initial separation";
  input Real M = 1.0 "Total mass";
  input Real eta = 0.25 "Symmetric mass ratio";
  output Real T "Coalescence time";
algorithm
  T := 5.0 / 256.0 * r0^4 / (eta * M^3);
  annotation(Documentation(info="<html>
<p><code>T = (5/256) r<sub>0</sub>&#8308;/(&eta;M&sup3;)</code>, the integral of
<a href=\"modelica://GeneralRelativity.CompactBinaries.inspiralRadialVelocity\">inspiralRadialVelocity</a>
(P. C. Peters, <i>Phys. Rev.</i> 136, B1224 (1964)). The time from <code>r<sub>0</sub></code> to <code>r</code> is
<code>(5/256)(r<sub>0</sub>&#8308; &minus; r&#8308;)/(&eta;M&sup3;)</code>.</p>
</html>"));
end mergerTimePeters;
