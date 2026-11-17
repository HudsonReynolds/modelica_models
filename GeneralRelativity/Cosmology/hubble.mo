within GeneralRelativity.Cosmology;
function hubble "Hubble rate H(a) in km/s/Mpc"
  input Real a "Scale factor (1 today)";
  input Parameters par;
  output Real H "Hubble rate, km/s/Mpc";
protected
  Real H0, Og, On, Ob, Oc, OL;
algorithm
  (H0, Og, On, Ob, Oc, OL) := densities(par);
  H := 100.0 * par.h * sqrt((Og + On) / a^4 + (Ob + Oc) / a^3 + par.omegaK / a^2 + OL);
  annotation(Documentation(info="<html>
<p>The Friedmann equation <code>H&sup2; = H<sub>0</sub>&sup2;(&Omega;<sub>r</sub>a<sup>&minus;4</sup> + &Omega;<sub>m</sub>a<sup>&minus;3</sup>
+ &Omega;<sub>k</sub>a<sup>&minus;2</sup> + &Omega;<sub>&Lambda;</sub>)</code>.</p>
</html>"));
end hubble;
