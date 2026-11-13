within GeneralRelativity.Connections;
function nonMetricity "Non-metricity Q_lmn = nabla_l g_mn of a connection with respect to a metric"
  input Real g[4, 4] "Covariant metric";
  input Real dg[4, 4, 4] "dg[l, m, n] = d_l g_mn";
  input Real Gamma[4, 4, 4] "Gamma[i, j, k] = Gamma^i_jk";
  output Real Q[4, 4, 4] "Q[l, m, n]; zero for a metric-compatible connection";
algorithm
  for l in 1:4 loop
    for m in 1:4 loop
      for n in 1:4 loop
        Q[l, m, n] := dg[l, m, n] - sum(Gamma[r, l, m] * g[r, n] + Gamma[r, l, n] * g[m, r] for r in 1:4);
      end for;
    end for;
  end for;
end nonMetricity;
