within GeneralRelativity;
package CompactBinaries "Compact binaries: post-Newtonian inspiral with radiation reaction, quadrupole waveform, merger remnant and ringdown"
  annotation(Documentation(info="<html>
<p>Two compact bodies (black holes or neutron stars) of masses <code>m<sub>1</sub>, m<sub>2</sub></code> in geometric
units <code>G = c = 1</code>. The relative orbit <code>x = x<sub>1</sub> &minus; x<sub>2</sub></code> moves under the
post-Newtonian relative acceleration in harmonic coordinates, with the 2.5PN radiation-reaction terms that
make the orbit shrink as gravitational waves carry energy away. Notation: total mass <code>M</code>, symmetric
mass ratio <code>&eta; = m<sub>1</sub>m<sub>2</sub>/M&sup2;</code>, <code>&gamma; = M/r</code>.</p>
<table>
<tr><th>Function</th><th>Gives</th></tr>
<tr><td><a href=\"modelica://GeneralRelativity.CompactBinaries.relativeAcceleration\">relativeAcceleration</a></td><td>Relative acceleration through 2PN (conservative) and 2.5PN (radiation reaction)</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.CompactBinaries.circularOmega\">circularOmega</a>, <a href=\"modelica://GeneralRelativity.CompactBinaries.circularOrbitSpeed\">circularOrbitSpeed</a>, <a href=\"modelica://GeneralRelativity.CompactBinaries.inspiralRadialVelocity\">inspiralRadialVelocity</a></td><td>Quasi-circular initial data: orbital frequency (Kepler's law with 2PN corrections) and the inward drift</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.CompactBinaries.mergerTimePeters\">mergerTimePeters</a></td><td>Time to coalescence of a circular orbit at leading order (Peters 1964)</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.CompactBinaries.quadrupoleStrain\">quadrupoleStrain</a></td><td>Gravitational-wave strain h<sub>+</sub>, h<sub>&times;</sub> for an observer on the orbital axis</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.CompactBinaries.finalSpin\">finalSpin</a>, <a href=\"modelica://GeneralRelativity.CompactBinaries.radiatedEnergyFraction\">radiatedEnergyFraction</a></td><td>Remnant black hole of non-spinning progenitors, from numerical-relativity fits</td></tr>
<tr><td><a href=\"modelica://GeneralRelativity.CompactBinaries.ringdownOmega\">ringdownOmega</a>, <a href=\"modelica://GeneralRelativity.CompactBinaries.ringdownQuality\">ringdownQuality</a></td><td>Fundamental l = m = 2 quasi-normal mode of the remnant: frequency and quality factor</td></tr>
</table>
<p><a href=\"modelica://GeneralRelativity.Examples.BinaryMerger\">Examples.BinaryMerger</a> puts these together:
inspiral, merger, remnant and ringdown, recorded for the viewer, whose ray tracer shows the two holes with
the superposed Kerr-Schild metric <a href=\"modelica://GeneralRelativity.Metrics.binaryKerrSchild\">Metrics.binaryKerrSchild</a>
and then the remnant Kerr hole. The post-Newtonian expansion is an expansion in <code>&gamma;</code> and
<code>v&sup2;</code>; it is accurate in the early inspiral and only qualitative in the last orbits before
merger (separations of a few <code>M</code>), where numerical relativity takes over.</p>
</html>"));
end CompactBinaries;
