within GeneralRelativity.Cosmology;
function scaleFactorAtTime "Scale factor at cosmic time t (Gyr), by Newton's method on cosmicTime"
  input Real t(unit = "Gyr");
  input Parameters par;
  output Real a;
protected
  Real H0, Og, On, Ob, Oc, OL;
  Real f;
  Real Hc;
algorithm
  (H0, Og, On, Ob, Oc, OL) := densities(par);
  // Start from the radiation-era solution a = (2 H0 sqrt(Or) t)^(1/2), capped.
  a := min(sqrt(2.0 * H0 * sqrt(Og + On) * t / MpcGyr), 1.0);
  for i in 1:60 loop
    f := cosmicTime(a, par) - t;
    Hc := a2H(a, H0, Og + On, Ob + Oc, par.omegaK, OL) / a "a H in 1/Mpc";
    // dt/da = 1/(a H)
    a := max(a - f / (MpcGyr / Hc), 0.5 * a);
    if abs(f) < 1.0e-12 * max(t, 1.0e-12) then
      break;
    end if;
  end for;
end scaleFactorAtTime;
