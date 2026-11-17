within GeneralRelativity.Cosmology;
function soundHorizon "Comoving sound horizon of the photon-baryon fluid, int_0^a c_s da/(a^2 H), in Mpc"
  input Real a;
  input Parameters par;
  input Integer n = 2000 "Simpson intervals";
  output Real rs(unit = "Mpc");
protected
  Real H0, Og, On, Ob, Oc, OL;
  Integer m = 2 * max(1, div(n, 2));
  Real dx, x, s;
algorithm
  (H0, Og, On, Ob, Oc, OL) := densities(par);
  // a = u^2, as in conformalTime.
  dx := sqrt(a) / m;
  s := 0.0;
  for i in 0:m loop
    x := i * dx;
    s := s + (if i == 0 or i == m then 1.0 else if mod(i, 2) == 1 then 4.0 else 2.0) * 2.0 * x
      / (a2H(x^2, H0, Og + On, Ob + Oc, par.omegaK, OL) * sqrt(3.0 * (1.0 + 0.75 * Ob / Og * x^2)));
  end for;
  rs := s * dx / 3.0;
  annotation(Documentation(info="<html>
<p>The sound speed of the tightly coupled photon-baryon fluid is <code>c<sub>s</sub> = c/&radic;(3(1 + R))</code>
with <code>R = 3&rho;<sub>b</sub>/(4&rho;<sub>&gamma;</sub>) &prop; a</code>. At last scattering
(<code>z<sub>*</sub> &asymp; 1090</code>) it is about 144 Mpc: the scale of the acoustic peaks.</p>
</html>"));
end soundHorizon;
