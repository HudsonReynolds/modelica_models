within GeneralRelativity.CompactBinaries;
function inspiralRadialVelocity "Inward drift dr/dt of a quasi-circular binary from gravitational-wave emission (leading order)"
  input Real r "Separation";
  input Real M = 1.0 "Total mass";
  input Real eta = 0.25 "Symmetric mass ratio";
  output Real rdot "dr/dt, negative";
algorithm
  rdot := -64.0 / 5.0 * eta * M^3 / r^3;
  annotation(Documentation(info="<html>
<p><code>dr/dt = &minus;(64/5)&eta;M&sup3;/r&sup3;</code> (Peters 1964): the energy <code>&minus;&eta;M&sup2;/(2r)</code> of a
circular orbit falls at the quadrupole rate. Starting an integration with this radial velocity, instead
of zero, avoids the residual eccentricity of a purely tangential start.</p>
</html>"));
end inspiralRadialVelocity;
