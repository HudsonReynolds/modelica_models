within GeneralRelativity.Cosmology;
function eventHorizon "Comoving distance light emitted at scale factor a can still travel, int_a^infinity da/(a^2 H), in Mpc"
  input Real a;
  input Parameters par;
  input Integer n = 2000 "Simpson intervals";
  output Real chi(unit = "Mpc");
protected
  Real H0, Og, On, Ob, Oc, OL;
  Integer m = 2 * max(1, div(n, 2));
  Real du, u, s;
algorithm
  (H0, Og, On, Ob, Oc, OL) := densities(par);
  assert(OL > 0.0, "eventHorizon: needs a cosmological constant (otherwise there is no event horizon)");
  // u = 1/a: int_0^{1/a} du / (H0 sqrt(Or u^4 + Om u^3 + Ok u^2 + OL)), smooth at u = 0.
  du := (1.0 / a) / m;
  s := 0.0;
  for i in 0:m loop
    u := i * du;
    s := s + (if i == 0 or i == m then 1.0 else if mod(i, 2) == 1 then 4.0 else 2.0) / (H0 * sqrt((Og + On) * u^4 + (Ob + Oc) * u^3 + par.omegaK * u^2 + OL));
  end for;
  chi := s * du / 3.0;
  annotation(Documentation(info="<html>
<p>Events beyond this comoving distance at scale factor <code>a</code> will never be seen by us: in a
universe whose expansion accelerates forever the integral converges. Multiply by <code>a</code> for the
proper distance at that time.</p>
</html>"));
end eventHorizon;
