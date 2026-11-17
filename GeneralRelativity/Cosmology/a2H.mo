within GeneralRelativity.Cosmology;
function a2H "da/d(eta) = a^2 H(a) in 1/Mpc (conformal time eta in Mpc); finite at a = 0"
  input Real a;
  input Real H0 "1/Mpc";
  input Real omegaR;
  input Real omegaM;
  input Real omegaK;
  input Real omegaLambda;
  output Real y;
algorithm
  y := H0 * sqrt(omegaR + omegaM * a + omegaK * a^2 + omegaLambda * a^4);
  annotation(Inline = true);
end a2H;
