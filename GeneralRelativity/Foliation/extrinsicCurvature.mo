within GeneralRelativity.Foliation;
function extrinsicCurvature "Extrinsic curvature K_ij of the slices x^0 = const from the metric and its first derivatives"
  input Real g[4, 4] "Covariant metric; x[1] = c*t";
  input Real dg[4, 4, 4] "dg[k, i, j] = d g_ij / d x^k";
  output Real K[3, 3] "K_ij = -(1/(2 alpha))(d_0 gamma_ij - D_i beta_j - D_j beta_i)";
protected
  Real alpha;
  Real betaUp[3];
  Real betaDown[3];
  Real gamma[3, 3];
  Real gammaInv[3, 3];
  Boolean accepted;
  Real Gamma3[3, 3, 3];
  Real Dbeta[3, 3] "Dbeta[i, j] = D_i beta_j";
algorithm
  (alpha, betaUp, betaDown, gamma, gammaInv, accepted) := admDecomposition(g);
  Gamma3 := GeneralRelativity.Tensors.christoffel(gammaInv, dg[2:4, 2:4, 2:4]);
  for i in 1:3 loop
    for j in 1:3 loop
      Dbeta[i, j] := dg[i + 1, 1, j + 1] - sum(Gamma3[k, i, j] * betaDown[k] for k in 1:3);
    end for;
  end for;
  for i in 1:3 loop
    for j in 1:3 loop
      K[i, j] := -(dg[1, i + 1, j + 1] - Dbeta[i, j] - Dbeta[j, i]) / (2.0 * alpha);
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>Uses the time derivative of the spatial metric and the spatial covariant
derivative of the lowered shift, both read from <code>dg</code>. The sign convention
(<code>K = -&nabla;<sub>a</sub>n<sup>a</sup></code>, negative for expanding
slices) is described in <a href=\"modelica://GeneralRelativity.Foliation\">Foliation</a>.
The result equals <code>-&alpha;&Gamma;<sup>0</sup><sub>ij</sub></code> of the
4-metric.</p>
</html>"));
end extrinsicCurvature;
