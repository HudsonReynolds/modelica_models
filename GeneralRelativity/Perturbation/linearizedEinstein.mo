within GeneralRelativity.Perturbation;
function linearizedEinstein "Linear-order Einstein tensor G_ab = R_ab - 1/2 eta_ab R from second derivatives of h"
  input Real ddh[4, 4, 4, 4] "ddh[k, l, a, b] = d^2 h_ab / (d x^k d x^l)";
  output Real G[4, 4] "G_ab to first order in h";
  output Real Ric[4, 4] "R_ab to first order";
  output Real ricciScalar "R to first order";
protected
  constant Real eta[4, 4] = [-1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
  Real R[4, 4, 4, 4];
algorithm
  R := linearizedRiemann(ddh);
  Ric := linearizedRicci(R);
  ricciScalar := sum(eta[a, a] * Ric[a, a] for a in 1:4);
  G := Ric - 0.5 * ricciScalar * eta;
  annotation(Documentation(info="<html>
<p>For the Newtonian-limit perturbation of a static potential this gives
<code>G<sub>00</sub> = 2&nabla;&sup2;&Phi;/c&sup2;</code> and <code>G<sub>ij</sub> = 0</code>, so
<code>G<sub>00</sub> = 8&pi;G&rho;/c&sup2;</code> is Poisson's equation.</p>
</html>"));
end linearizedEinstein;
