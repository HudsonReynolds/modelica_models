within GeneralRelativity.Cosmology;
function cosmicTime "Cosmic (proper) time since the big bang at scale factor a, in Gyr"
  input Real a "Scale factor";
  input Parameters par;
  input Integer n = 2000 "Simpson intervals";
  output Real t(unit = "Gyr");
protected
  Real H0, Og, On, Ob, Oc, OL;
  Integer m = 2 * max(1, div(n, 2));
  Real dx, x, s;
algorithm
  (H0, Og, On, Ob, Oc, OL) := densities(par);
  // dt = da/(a H) = a da/(a^2 H); with a = u^2, dt = 2u^3 du/(a^2 H), smooth at u = 0.
  dx := sqrt(a) / m;
  s := 0.0;
  for i in 0:m loop
    x := i * dx;
    s := s + (if i == 0 or i == m then 1.0 else if mod(i, 2) == 1 then 4.0 else 2.0) * 2.0 * x^3 / a2H(x^2, H0, Og + On, Ob + Oc, par.omegaK, OL);
  end for;
  t := s * dx / 3.0 * MpcGyr;
  annotation(Documentation(info="<html>
<p><code>t(a) = &int;<sub>0</sub><sup>a</sup> da'/(a'H)</code>, Simpson's rule in <code>u = &radic;a</code> (the integrand
<code>2u&sup3;/(a&sup2;H)</code> is smooth from the big bang on and resolves the radiation-matter transition). With the Planck 2018 parameters the age is 13.81 Gyr
(13.80 with Planck's massive neutrino).</p>
</html>"));
end cosmicTime;
