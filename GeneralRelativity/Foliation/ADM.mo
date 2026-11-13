within GeneralRelativity.Foliation;
record ADM "3+1 split of a metric at one event: lapse, shift, spatial metric, extrinsic curvature, 3-curvature and constraint residuals"
  Real x[4] "Coordinates of the event; x[1] = c*t";
  Real alpha "Lapse alpha";
  Real betaUp[3] "Shift beta^i";
  Real betaDown[3] "Shift beta_i = g_0i";
  Real gamma[3, 3] "Spatial metric gamma_ij = g_ij";
  Real gammaInv[3, 3] "Inverse spatial metric gamma^ij";
  Real nUp[4] "Future unit normal n^a = (1, -beta^i)/alpha";
  Real nDown[4] "Unit normal n_a = (-alpha, 0, 0, 0)";
  Real K[3, 3] "Extrinsic curvature K_ij = -(1/(2 alpha))(d_0 gamma_ij - D_i beta_j - D_j beta_i), d_0 = d/d(ct)";
  Real KMixed[3, 3] "Mixed extrinsic curvature K^i_j = gamma^ik K_kj";
  Real trK "Mean curvature K = gamma^ij K_ij = -div n; negative for expanding slices";
  Real Gamma3[3, 3, 3] "Spatial Christoffel symbols Gamma3[i, j, k] = (3)Gamma^i_jk";
  Real ricci3[3, 3] "Spatial Ricci tensor (3)R_ij";
  Real ricciScalar3 "Spatial Ricci scalar (3)R";
  Real normalEinstein "G_ab n^a n^b = 8 pi G rho / c^4 for Eulerian observers";
  Real momentumEinstein[3] "-G_ab n^a delta^b_i = 8 pi G j_i / c^4 for Eulerian observers";
  Real hamiltonian "Gauss residual (3)R + K^2 - K_ij K^ij - 2 G_ab n^a n^b; zero for any metric";
  Real momentum[3] "Codazzi residual D_j K^j_i - D_i K + G_ab n^a delta^b_i; zero for any metric";
  Real valid "1 when the slice is spacelike and both metrics were inverted, else 0";
  annotation(Documentation(info="<html>
<p>Filled by <a href=\"modelica://GeneralRelativity.Foliation.admAt\">admAt</a> or
<a href=\"modelica://GeneralRelativity.Foliation.admFromJet\">admFromJet</a>.
Conventions and the meaning of the residuals are in
<a href=\"modelica://GeneralRelativity.Foliation\">Foliation</a>. Spatial index
<code>i = 1..3</code> is spacetime index <code>i + 1</code>.</p>
</html>"));
end ADM;
