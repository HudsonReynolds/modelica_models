within GeneralRelativity.CompactBinaries;
function circularOmega "Orbital angular frequency of a circular binary at harmonic separation r, Kepler's law with 1PN and 2PN corrections"
  input Real r "Harmonic-coordinate separation";
  input Real M = 1.0 "Total mass";
  input Real eta = 0.25 "Symmetric mass ratio";
  input Integer pnOrder = 2 "0, 1 or 2";
  output Real omega "Orbital angular frequency d phi / dt";
protected
  Real gamma = M / r;
  Real w2 = 1.0;
algorithm
  if pnOrder >= 1 then
    w2 := w2 + (-3.0 + eta) * gamma;
  end if;
  if pnOrder >= 2 then
    w2 := w2 + (6.0 + 41.0 / 4.0 * eta + eta^2) * gamma^2;
  end if;
  omega := sqrt(M / r^3 * w2);
  annotation(Documentation(info="<html>
<p><code>&Omega;&sup2; = (M/r&sup3;)[1 + (&minus;3 + &eta;)&gamma; + (6 + 41&eta;/4 + &eta;&sup2;)&gamma;&sup2;]</code>, harmonic
coordinates (Blanchet, <i>Living Rev. Relativ.</i> 17, 2 (2014)). It is the conservative part of
<a href=\"modelica://GeneralRelativity.CompactBinaries.relativeAcceleration\">relativeAcceleration</a> on a
circular orbit.</p>
</html>"));
end circularOmega;
