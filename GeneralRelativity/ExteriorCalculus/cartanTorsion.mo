within GeneralRelativity.ExteriorCalculus;
function cartanTorsion "Residual of Cartan's first structure equation, T^a = d theta^a + omega^a_b ^ theta^b (zero for Levi-Civita)"
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric function";
  input TetradField tetrad "Orthonormal tetrad field for this metric";
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Metric parameters";
  input Real q[:] "Tetrad parameters";
  input Real step[4] "Finite-difference step per coordinate";
  output Real T[4, 4, 4] "T[a, mu, nu] = d_mu theta^a_nu - d_nu theta^a_mu + omega^a_b mu theta^b_nu - omega^a_b nu theta^b_mu";
protected
  Real omega[4, 4, 4];
  Real theta[4, 4];
  Real thetaPlus[4, 4];
  Real thetaMinus[4, 4];
  Real dtheta[4, 4, 4] "dtheta[mu, a, nu] = d_mu theta^a_nu";
  Real xPlus[4];
  Real xMinus[4];
  Real gPlus[4, 4];
  Real gMinus[4, 4];
  Real ePlus[4, 4];
  Real eMinus[4, 4];
  Real e[4, 4];
  Real g[4, 4];
algorithm
  omega := spinConnection(metric, tetrad, x, p, q, step);
  g := metric(x, p);
  e := tetrad(x, q);
  theta := coframe(g, e);
  for mu in 1:4 loop
    xPlus := x;
    xMinus := x;
    xPlus[mu] := x[mu] + step[mu];
    xMinus[mu] := x[mu] - step[mu];
    gPlus := metric(xPlus, p);
    gMinus := metric(xMinus, p);
    ePlus := tetrad(xPlus, q);
    eMinus := tetrad(xMinus, q);
    thetaPlus := coframe(gPlus, ePlus);
    thetaMinus := coframe(gMinus, eMinus);
    dtheta[mu, :, :] := (thetaPlus - thetaMinus) / (2.0 * step[mu]);
  end for;
  for a in 1:4 loop
    for mu in 1:4 loop
      for nu in 1:4 loop
        T[a, mu, nu] := dtheta[mu, a, nu] - dtheta[nu, a, mu]
          + sum(omega[a, b, mu] * theta[b, nu] - omega[a, b, nu] * theta[b, mu] for b in 1:4);
      end for;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>For the Levi-Civita connection the torsion vanishes, so the result is a
finite-difference residual of order <code>step&sup2;</code>. A nonzero value
signals a tetrad that is not orthonormal for <code>metric</code> or a wrong
connection.</p>
</html>"));
end cartanTorsion;
