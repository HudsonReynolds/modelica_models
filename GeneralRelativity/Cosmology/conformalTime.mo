within GeneralRelativity.Cosmology;
function conformalTime "Conformal time eta(a) = int_0^a da'/(a'^2 H) in Mpc: the comoving distance light travels since the big bang"
  input Real a "Scale factor";
  input Parameters par;
  input Integer n = 2000 "Simpson intervals";
  output Real eta(unit = "Mpc");
protected
  Real H0, Og, On, Ob, Oc, OL;
  Integer m = 2 * max(1, div(n, 2));
  Real dx, x, s;
algorithm
  (H0, Og, On, Ob, Oc, OL) := densities(par);
  // a = u^2: d(eta) = 2u du / (a^2 H(a)), smooth from u = 0 and resolving the radiation-matter transition.
  dx := sqrt(a) / m;
  s := 0.0;
  for i in 0:m loop
    x := i * dx;
    s := s + (if i == 0 or i == m then 1.0 else if mod(i, 2) == 1 then 4.0 else 2.0) * 2.0 * x / a2H(x^2, H0, Og + On, Ob + Oc, par.omegaK, OL);
  end for;
  eta := s * dx / 3.0;
end conformalTime;
