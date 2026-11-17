within GeneralRelativity.CompactBinaries;
function radiatedEnergyFraction "Fraction of the total mass radiated in gravitational waves by the merger of two non-spinning black holes"
  input Real eta "Symmetric mass ratio";
  output Real e "E_rad / M";
protected
  constant Real p0 = 0.04827;
  constant Real eIsco = sqrt(8.0 / 9.0) "Specific energy at the Schwarzschild innermost stable circular orbit";
algorithm
  e := (1.0 - eIsco) * eta + 4.0 * eta^2 * (4.0 * p0 + eIsco - 1.0);
  annotation(Documentation(info="<html>
<p>E. Barausse, V. Morozova and L. Rezzolla, <i>Astrophys. J.</i> 758, 63 (2012), with zero progenitor
spin: <code>E/M = (1 &minus; &Etilde;<sub>ISCO</sub>)&eta; + 4&eta;&sup2;(4p<sub>0</sub> + &Etilde;<sub>ISCO</sub> &minus; 1)</code>,
<code>&Etilde;<sub>ISCO</sub> = &radic;(8/9)</code>, <code>p<sub>0</sub> = 0.04827</code>. Equal masses radiate 4.8% of
their mass; the remnant mass is <code>M(1 &minus; E/M)</code>.</p>
</html>"));
end radiatedEnergyFraction;
