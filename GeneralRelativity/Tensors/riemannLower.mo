within GeneralRelativity.Tensors;
function riemannLower "Fully covariant Riemann tensor R[a, b, c, d] = R_abcd from metric second derivatives and Christoffel symbols"
  input Real g[4, 4] "Covariant metric";
  input Real ddg[4, 4, 4, 4] "ddg[k, l, i, j] = d^2 g_ij / (d x^k d x^l)";
  input Real Gamma[4, 4, 4] "Gamma[i, j, k] = Gamma^i_jk";
  output Real R[4, 4, 4, 4] "R[a, b, c, d] = R_abcd";
protected
  Real quadratic;
algorithm
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        for d in 1:4 loop
          quadratic := 0.0;
          for e in 1:4 loop
            for f in 1:4 loop
              quadratic := quadratic + g[e, f] * (Gamma[e, b, c] * Gamma[f, a, d] - Gamma[e, b, d] * Gamma[f, a, c]);
            end for;
          end for;
          R[a, b, c, d] := 0.5 * (ddg[b, c, a, d] + ddg[a, d, b, c] - ddg[a, c, b, d] - ddg[b, d, a, c]) + quadratic;
        end for;
      end for;
    end for;
  end for;
end riemannLower;
