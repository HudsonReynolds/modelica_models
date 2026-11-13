within GeneralRelativity.Foliation;
function admFromJet "ADM record from a metric and its first and second derivatives at one event"
  input Real x[4] "Coordinates of the event; stored in the record";
  input Real g[4, 4] "Covariant metric";
  input Real dg[4, 4, 4] "dg[k, i, j] = d g_ij / d x^k";
  input Real ddg[4, 4, 4, 4] "ddg[k, l, i, j] = d^2 g_ij / (d x^k d x^l)";
  output ADM s "3+1 quantities and constraint residuals";
protected
  Boolean splitAccepted;
  Boolean accepted4;
  Real gInv[4, 4];
  Real Gamma4[4, 4, 4];
  Real riemann4[4, 4, 4, 4];
  Real ricci4[4, 4];
  Real einstein4[4, 4];
  Real riemann3[3, 3, 3, 3];
  Real KUp[3, 3] "K^ij";
  Real dK[4, 3, 3];
  Real DK[3, 3, 3] "DK[k, i, j] = D_k K_ij";
  Real divergence;
  Real gradient;
algorithm
  s.x := x;
  (s.alpha, s.betaUp, s.betaDown, s.gamma, s.gammaInv, splitAccepted) := admDecomposition(g);
  (s.nUp, s.nDown) := normalVector(g);
  (gInv, accepted4) := GeneralRelativity.Tensors.inverseMetric(g);
  s.valid := if splitAccepted and accepted4 then 1.0 else 0.0;

  // Extrinsic curvature and its spatial covariant derivative.
  s.K := extrinsicCurvature(g, dg);
  s.KMixed := s.gammaInv * s.K;
  s.trK := sum(s.KMixed[i, i] for i in 1:3);
  KUp := s.KMixed * s.gammaInv;

  // Intrinsic curvature of the slice from the spatial block of the jet.
  s.Gamma3 := GeneralRelativity.Tensors.christoffel(s.gammaInv, dg[2:4, 2:4, 2:4]);
  riemann3 := GeneralRelativity.Tensors.riemannLower(s.gamma, ddg[2:4, 2:4, 2:4, 2:4], s.Gamma3);
  s.ricci3 := GeneralRelativity.Tensors.ricci(s.gammaInv, riemann3);
  s.ricciScalar3 := GeneralRelativity.Tensors.ricciScalar(s.gammaInv, s.ricci3);

  // Spacetime Einstein tensor from the same jet.
  Gamma4 := GeneralRelativity.Tensors.christoffel(gInv, dg);
  riemann4 := GeneralRelativity.Tensors.riemannLower(g, ddg, Gamma4);
  ricci4 := GeneralRelativity.Tensors.ricci(gInv, riemann4);
  einstein4 := GeneralRelativity.Tensors.einstein(g, ricci4, GeneralRelativity.Tensors.ricciScalar(gInv, ricci4));
  s.normalEinstein := s.nUp * einstein4 * s.nUp;
  for i in 1:3 loop
    s.momentumEinstein[i] := -sum(s.nUp[a] * einstein4[a, i + 1] for a in 1:4);
  end for;

  // Gauss: (3)R + K^2 - K_ij K^ij = 2 G_ab n^a n^b.
  s.hamiltonian := s.ricciScalar3 + s.trK^2 - sum(s.K[i, j] * KUp[i, j] for i in 1:3, j in 1:3) - 2.0 * s.normalEinstein;

  // Codazzi: D_j K^j_i - D_i K = -G_ab n^a delta^b_i.
  dK := extrinsicCurvatureGradient(g, dg, ddg);
  for k in 1:3 loop
    for i in 1:3 loop
      for j in 1:3 loop
        DK[k, i, j] := dK[k + 1, i, j] - sum(s.Gamma3[l, k, i] * s.K[l, j] + s.Gamma3[l, k, j] * s.K[i, l] for l in 1:3);
      end for;
    end for;
  end for;
  for i in 1:3 loop
    divergence := sum(s.gammaInv[j, k] * DK[k, j, i] for j in 1:3, k in 1:3);
    gradient := sum(s.gammaInv[j, k] * DK[i, j, k] for j in 1:3, k in 1:3);
    s.momentum[i] := divergence - gradient - s.momentumEinstein[i];
  end for;
  annotation(Documentation(info="<html>
<p>Pure algebra on the jet <code>(g, &part;g, &part;&part;g)</code>: analytic
derivatives work as well as finite differences. The 3-curvature uses the
dimension-generic <a href=\"modelica://GeneralRelativity.Tensors\">Tensors</a>
chain on the spatial block; the 4-dimensional Einstein tensor uses the full
jet. Because the Gauss and Codazzi equations are identities on the 2-jet, the
returned residuals are rounding error for any input; their size relative to
the individual terms measures the numerical health of the split. See
<a href=\"modelica://GeneralRelativity.Foliation\">Foliation</a> for the formulas.</p>
</html>"));
end admFromJet;
