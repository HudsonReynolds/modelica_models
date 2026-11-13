within GeneralRelativity.ExteriorCalculus;
function spinConnection "Levi-Civita spin connection omega[a, b, mu] = omega^a_b mu of a tetrad field"
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric function";
  input TetradField tetrad "Orthonormal tetrad field for this metric";
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Metric parameters";
  input Real q[:] "Tetrad parameters";
  input Real step[4] "Finite-difference step per coordinate";
  output Real omega[4, 4, 4] "omega^a_b mu = theta^a_nu (d_mu e_b^nu + Gamma^nu_mu lambda e_b^lambda)";
protected
  Real g[4, 4];
  Real gInv[4, 4];
  Boolean accepted;
  Real dg[4, 4, 4];
  Real Gamma[4, 4, 4];
  Real e[4, 4];
  Real theta[4, 4];
  Real de[4, 4, 4] "de[mu, b, nu] = d_mu e_b^nu";
  Real xPlus[4];
  Real xMinus[4];
  Real ePlus[4, 4];
  Real eMinus[4, 4];
  Real covariant;
algorithm
  g := metric(x, p);
  (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(g);
  dg := GeneralRelativity.Metrics.derivative(metric, x, p, step);
  Gamma := GeneralRelativity.Tensors.christoffel(gInv, dg);
  e := tetrad(x, q);
  theta := coframe(g, e);
  for mu in 1:4 loop
    xPlus := x;
    xMinus := x;
    xPlus[mu] := x[mu] + step[mu];
    xMinus[mu] := x[mu] - step[mu];
    ePlus := tetrad(xPlus, q);
    eMinus := tetrad(xMinus, q);
    de[mu, :, :] := (ePlus - eMinus) / (2.0 * step[mu]);
  end for;
  for a in 1:4 loop
    for b in 1:4 loop
      for mu in 1:4 loop
        omega[a, b, mu] := 0.0;
        for nu in 1:4 loop
          covariant := de[mu, b, nu] + sum(Gamma[nu, mu, l] * e[b, l] for l in 1:4);
          omega[a, b, mu] := omega[a, b, mu] + theta[a, nu] * covariant;
        end for;
      end for;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>Defined by <code>&nabla;<sub>&mu;</sub>e<sub>b</sub> = &omega;<sup>a</sup><sub>b&mu;</sub>e<sub>a</sub></code>,
so <code>omega[a, b, :]</code> is the connection 1-form
<code>&omega;<sup>a</sup><sub>b</sub></code> in coordinate components. Metric
compatibility makes <code>&omega;<sub>ab</sub> = &eta;<sub>ac</sub>&omega;<sup>c</sup><sub>b</sub></code>
antisymmetric, and the first structure equation
<code>d&theta;<sup>a</sup> + &omega;<sup>a</sup><sub>b</sub>&and;&theta;<sup>b</sup> = 0</code>
holds (<a href=\"modelica://GeneralRelativity.ExteriorCalculus.cartanTorsion\">cartanTorsion</a>).
Metric and tetrad derivatives are central differences with <code>step</code>.</p>
</html>"));
end spinConnection;
