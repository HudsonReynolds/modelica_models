within GeneralRelativity.Tensors;
record Curvature "Metric and curvature tensors at one event"
  Real x[4] "Coordinates of the event; x[1] = c*t";
  Real g[4, 4] "Covariant metric g_ab";
  Real gInv[4, 4] "Contravariant metric g^ab";
  Real Gamma[4, 4, 4] "Christoffel symbols Gamma[i, j, k] = Gamma^i_jk";
  Real riemann[4, 4, 4, 4] "R_abcd";
  Real ricci[4, 4] "R_ab";
  Real ricciScalar "R";
  Real einstein[4, 4] "G_ab";
  Real kretschmann "R_abcd R^abcd";
  Real valid "1 when the metric was inverted at x, 0 when it is numerically singular there; a Real so the record can sit in continuous equations";
  annotation(Documentation(info="<html>
<p>Filled by <a href=\"modelica://GeneralRelativity.Metrics.curvatureAt\">Metrics.curvatureAt</a>.
The implied stress-energy follows from
<a href=\"modelica://GeneralRelativity.Tensors.stressEnergyFromEinstein\">stressEnergyFromEinstein</a>
once <code>G</code> and <code>c</code> are chosen.</p>
</html>"));
end Curvature;
