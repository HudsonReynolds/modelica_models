within GeneralRelativity.Cosmology;
function decelerationParameter "Deceleration parameter q = -a a''/a'^2 at scale factor a"
  input Real a;
  input Parameters par;
  output Real q;
protected
  Real H0, Og, On, Ob, Oc, OL;
  Real E2;
algorithm
  (H0, Og, On, Ob, Oc, OL) := densities(par);
  E2 := (Og + On) / a^4 + (Ob + Oc) / a^3 + par.omegaK / a^2 + OL;
  // q = sum_i Omega_i(a) (1 + 3 w_i) / 2
  q := ((Og + On) / a^4 + 0.5 * (Ob + Oc) / a^3 - OL) / E2;
end decelerationParameter;
