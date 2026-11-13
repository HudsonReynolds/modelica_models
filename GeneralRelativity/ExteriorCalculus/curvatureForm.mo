within GeneralRelativity.ExteriorCalculus;
function curvatureForm "Curvature 2-forms from Cartan's second structure equation, Omega^a_b = d omega^a_b + omega^a_c ^ omega^c_b"
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric function";
  input TetradField tetrad "Orthonormal tetrad field for this metric";
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Metric parameters";
  input Real q[:] "Tetrad parameters";
  input Real step[4] "Outer finite-difference step; the spin connection uses step/2";
  output Real Omega[4, 4, 4, 4] "Omega[a, b, mu, nu] = R^a_b mu nu";
protected
  Real omega[4, 4, 4];
  Real omegaPlus[4, 4, 4];
  Real omegaMinus[4, 4, 4];
  Real domega[4, 4, 4, 4] "domega[mu, a, b, nu] = d_mu omega^a_b nu";
  Real xPlus[4];
  Real xMinus[4];
  Real innerStep[4] = 0.5 * step;
algorithm
  omega := spinConnection(metric, tetrad, x, p, q, innerStep);
  for mu in 1:4 loop
    xPlus := x;
    xMinus := x;
    xPlus[mu] := x[mu] + step[mu];
    xMinus[mu] := x[mu] - step[mu];
    omegaPlus := spinConnection(metric, tetrad, xPlus, p, q, innerStep);
    omegaMinus := spinConnection(metric, tetrad, xMinus, p, q, innerStep);
    for a in 1:4 loop
      for b in 1:4 loop
        for nu in 1:4 loop
          domega[mu, a, b, nu] := (omegaPlus[a, b, nu] - omegaMinus[a, b, nu]) / (2.0 * step[mu]);
        end for;
      end for;
    end for;
  end for;
  for a in 1:4 loop
    for b in 1:4 loop
      for mu in 1:4 loop
        for nu in 1:4 loop
          Omega[a, b, mu, nu] := domega[mu, a, b, nu] - domega[nu, a, b, mu]
            + sum(omega[a, c, mu] * omega[c, b, nu] - omega[a, c, nu] * omega[c, b, mu] for c in 1:4);
        end for;
      end for;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>The exterior derivative of the connection forms differences
<a href=\"modelica://GeneralRelativity.ExteriorCalculus.spinConnection\">spinConnection</a>
at neighbouring events (nested central differences: 9 connection evaluations of
9 metric and 9 tetrad evaluations each). The inner step is half the outer one so
that the inner rounding error, of order <code>&epsilon;/step</code>, is not
amplified by the outer division more than necessary while the stencils stay
distinct; the total error is of order <code>step&sup2;</code>. With a
correct connection, <code>Omega</code> equals
<a href=\"modelica://GeneralRelativity.ExteriorCalculus.frameRiemann\">frameRiemann</a>
of the coordinate Riemann tensor.</p>
</html>"));
end curvatureForm;
