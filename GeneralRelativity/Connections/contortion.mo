within GeneralRelativity.Connections;
function contortion "Contortion K^i_jk of a metric-compatible connection with torsion T: Gamma = Levi-Civita + K"
  input Real g[4, 4] "Covariant metric";
  input Real gInv[4, 4] "Contravariant metric";
  input Real T[4, 4, 4] "Torsion T^i_jk";
  output Real K[4, 4, 4] "K^i_jk";
protected
  Real TLower[4, 4, 4] "T_ajk = g_ai T^i_jk";
  Real KLower[4, 4, 4] "K_ajk";
algorithm
  for a in 1:4 loop
    for j in 1:4 loop
      for k in 1:4 loop
        TLower[a, j, k] := sum(g[a, i] * T[i, j, k] for i in 1:4);
      end for;
    end for;
  end for;
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        KLower[a, b, c] := 0.5 * (TLower[a, b, c] - TLower[b, c, a] + TLower[c, a, b]);
      end for;
    end for;
  end for;
  for i in 1:4 loop
    for j in 1:4 loop
      for k in 1:4 loop
        K[i, j, k] := sum(gInv[i, a] * KLower[a, j, k] for a in 1:4);
      end for;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>Metric compatibility makes <code>K<sub>abc</sub></code> antisymmetric in its first and
last index, and <code>T<sub>abc</sub> = K<sub>abc</sub> &minus; K<sub>acb</sub></code>; the unique
solution is <code>K<sub>abc</sub> = &frac12;(T<sub>abc</sub> &minus; T<sub>bca</sub> + T<sub>cab</sub>)</code>
with all indices in the order of <code>&Gamma;<sup>a</sup><sub>bc</sub></code>. Adding <code>K</code> to the
Christoffel symbols gives the metric connection with torsion <code>T</code>
(Einstein-Cartan, teleparallel).</p>
</html>"));
end contortion;
