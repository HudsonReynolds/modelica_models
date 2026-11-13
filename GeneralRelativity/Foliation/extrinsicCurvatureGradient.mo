within GeneralRelativity.Foliation;
function extrinsicCurvatureGradient "Coordinate gradient dK[k, i, j] = d K_ij / d x^k from the metric 2-jet, without nested differencing"
  input Real g[4, 4] "Covariant metric; x[1] = c*t";
  input Real dg[4, 4, 4] "dg[k, i, j] = d g_ij / d x^k";
  input Real ddg[4, 4, 4, 4] "ddg[k, l, i, j] = d^2 g_ij / (d x^k d x^l)";
  output Real dK[4, 3, 3] "dK[k, i, j] = d_k K_ij; k is a spacetime index, i and j spatial";
protected
  Real gInv[4, 4];
  Boolean accepted;
  Real alpha;
  Real nUp[4];
  Real lower[4, 3, 3] "lower[l, i, j] = Gamma_l(i+1)(j+1), Christoffel symbols of the first kind";
  Real dgInv[4, 4, 4] "dgInv[k, a, b] = d_k g^ab";
  Real dAlpha[4];
  Real dn[4, 4] "dn[k, l] = d_k n^l";
  Real dLower;
algorithm
  (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(g);
  alpha := 1.0 / sqrt(max(-gInv[1, 1], Modelica.Constants.small));
  nUp := -alpha * gInv[:, 1];
  for l in 1:4 loop
    for i in 1:3 loop
      for j in 1:3 loop
        lower[l, i, j] := 0.5 * (dg[i + 1, l, j + 1] + dg[j + 1, l, i + 1] - dg[l, i + 1, j + 1]);
      end for;
    end for;
  end for;
  for k in 1:4 loop
    dgInv[k, :, :] := -gInv * dg[k, :, :] * gInv;
    dAlpha[k] := 0.5 * alpha^3 * dgInv[k, 1, 1];
    for l in 1:4 loop
      dn[k, l] := -dAlpha[k] * gInv[l, 1] - alpha * dgInv[k, l, 1];
    end for;
  end for;
  for k in 1:4 loop
    for i in 1:3 loop
      for j in 1:3 loop
        dK[k, i, j] := 0.0;
        for l in 1:4 loop
          dLower := 0.5 * (ddg[k, i + 1, l, j + 1] + ddg[k, j + 1, l, i + 1] - ddg[k, l, i + 1, j + 1]);
          dK[k, i, j] := dK[k, i, j] + dn[k, l] * lower[l, i, j] + nUp[l] * dLower;
        end for;
      end for;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>Differentiates <code>K<sub>ij</sub> = n<sup>l</sup>&Gamma;<sub>lij</sub></code> with
<code>n<sup>l</sup> = -&alpha;g<sup>l0</sup></code>,
<code>&alpha; = (-g<sup>00</sup>)<sup>-1/2</sup></code> and
<code>&part;<sub>k</sub>g<sup>ab</sup> = -g<sup>ac</sup>&part;<sub>k</sub>g<sub>cd</sub>g<sup>db</sup></code>
by the chain rule. Only the metric and its first and second derivatives at one
event are needed, so the momentum constraint costs no extra metric evaluations
and, being an algebraic identity in the jet, holds to rounding error even with
finite-difference derivatives. Differencing
<a href=\"modelica://GeneralRelativity.Foliation.extrinsicCurvature\">extrinsicCurvature</a>
at neighbouring events gives the same result up to truncation error.</p>
</html>"));
end extrinsicCurvatureGradient;
