within GeneralRelativity.Cosmology;
function comovingDistance "Comoving distance travelled by light between scale factors a1 < a2, in Mpc"
  input Real a1;
  input Real a2;
  input Parameters par;
  input Integer n = 2000 "Simpson intervals";
  output Real chi(unit = "Mpc");
protected
  Real H0, Og, On, Ob, Oc, OL;
  Integer m = 2 * max(1, div(n, 2));
  Real dx, x, s;
algorithm
  (H0, Og, On, Ob, Oc, OL) := densities(par);
  // a = u^2, as in conformalTime.
  dx := (sqrt(a2) - sqrt(a1)) / m;
  s := 0.0;
  for i in 0:m loop
    x := sqrt(a1) + i * dx;
    s := s + (if i == 0 or i == m then 1.0 else if mod(i, 2) == 1 then 4.0 else 2.0) * 2.0 * x / a2H(x^2, H0, Og + On, Ob + Oc, par.omegaK, OL);
  end for;
  chi := s * dx / 3.0;
end comovingDistance;
