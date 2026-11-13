within GeneralRelativity.Tensors;
function riemannLower "Fully covariant Riemann tensor R[a, b, c, d] = R_abcd from metric second derivatives and Christoffel symbols, in any dimension"
  input Real g[:, size(g, 1)] "Covariant metric";
  input Real ddg[size(g, 1), size(g, 1), size(g, 1), size(g, 1)] "ddg[k, l, i, j] = d^2 g_ij / (d x^k d x^l)";
  input Real Gamma[size(g, 1), size(g, 1), size(g, 1)] "Gamma[i, j, k] = Gamma^i_jk";
  output Real R[size(g, 1), size(g, 1), size(g, 1), size(g, 1)] "R[a, b, c, d] = R_abcd";
protected
  Integer n = size(g, 1);
  Real quadratic;
algorithm
  for a in 1:n loop
    for b in 1:n loop
      for c in 1:n loop
        for d in 1:n loop
          quadratic := 0.0;
          for e in 1:n loop
            for f in 1:n loop
              quadratic := quadratic + g[e, f] * (Gamma[e, b, c] * Gamma[f, a, d] - Gamma[e, b, d] * Gamma[f, a, c]);
            end for;
          end for;
          R[a, b, c, d] := 0.5 * (ddg[b, c, a, d] + ddg[a, d, b, c] - ddg[a, c, b, d] - ddg[b, d, a, c]) + quadratic;
        end for;
      end for;
    end for;
  end for;
end riemannLower;
