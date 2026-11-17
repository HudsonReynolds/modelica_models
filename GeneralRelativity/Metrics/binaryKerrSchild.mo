within GeneralRelativity.Metrics;
function binaryKerrSchild "Two Schwarzschild black holes superposed in Kerr-Schild form, Cartesian {ct, x, y, z}; p = {m1, m2, x1, y1, z1, x2, y2, z2}"
  extends GeneralRelativity.Interfaces.MetricFunction;
protected
  Real r1 = sqrt((x[2] - p[3])^2 + (x[3] - p[4])^2 + (x[4] - p[5])^2) "Distance to body 1";
  Real r2 = sqrt((x[2] - p[6])^2 + (x[3] - p[7])^2 + (x[4] - p[8])^2) "Distance to body 2";
  Real k1[4] = {1.0, (x[2] - p[3]) / r1, (x[3] - p[4]) / r1, (x[4] - p[5]) / r1};
  Real k2[4] = {1.0, (x[2] - p[6]) / r2, (x[3] - p[7]) / r2, (x[4] - p[8]) / r2};
  Real f1 = 2.0 * p[1] / r1;
  Real f2 = 2.0 * p[2] / r2;
algorithm
  g := diagonal({-1.0, 1.0, 1.0, 1.0});
  for i in 1:4 loop
    for j in 1:4 loop
      g[i, j] := g[i, j] + f1 * k1[i] * k1[j] + f2 * k2[i] * k2[j];
    end for;
  end for;
  annotation(Documentation(info="<html>
<p><code>g<sub>&mu;&nu;</sub> = &eta;<sub>&mu;&nu;</sub> + &Sigma;<sub>a</sub> (2m<sub>a</sub>/r<sub>a</sub>) k<sup>(a)</sup><sub>&mu;</sub>k<sup>(a)</sup><sub>&nu;</sub></code>,
with <code>k<sup>(a)</sup> = (1, (x &minus; x<sub>a</sub>)/r<sub>a</sub>)</code> and <code>r<sub>a</sub> = |x &minus; x<sub>a</sub>|</code>:
the Kerr-Schild form of <a href=\"modelica://GeneralRelativity.Metrics.kerrSchild\">kerrSchild</a> (with
<code>a = 0</code>) written once for each body and added. Positions are the bodies' positions at the
moment of the slice.</p>
<p>This <b>superposed Kerr-Schild</b> metric is not a solution of the field equations: it is exact for
one body (either mass zero) and far from both, and is the classic starting point for binary initial data
(R. A. Matzner, M. F. Huq and D. Shoemaker, <i>Phys. Rev. D</i> 59, 024015 (1998)). Each term is regular across
its horizon, so null geodesics can be integrated straight into either hole. The viewer's ray tracer
uses exactly this formula for the <code>binary</code> metric family (with the bodies at their current
positions, frozen while a ray is traced), switching to the remnant Kerr hole after the merger; see
<a href=\"modelica://GeneralRelativity.Examples.BinaryMerger\">Examples.BinaryMerger</a>.</p>
</html>"));
end binaryKerrSchild;
