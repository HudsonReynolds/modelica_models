within GeneralRelativity.Cosmology;
function densities "Density parameters today and the Hubble constant in units of 1/Mpc"
  input Parameters par;
  output Real H0(unit = "1/Mpc") "Hubble constant divided by c";
  output Real omegaGamma "Photons";
  output Real omegaNu "Massless neutrinos";
  output Real omegaB "Baryons";
  output Real omegaC "Cold dark matter";
  output Real omegaLambda "Cosmological constant, from closure";
algorithm
  H0 := par.h * 100.0 / (GeneralRelativity.Constants.c / 1000.0);
  omegaGamma := 2.47282e-5 * (par.TCMB / 2.7255)^4 / par.h^2;
  omegaNu := omegaGamma * par.Neff * 7.0 / 8.0 * (4.0 / 11.0)^(4.0 / 3.0);
  omegaB := par.omegaBh2 / par.h^2;
  omegaC := par.omegaCh2 / par.h^2;
  omegaLambda := 1.0 - par.omegaK - omegaGamma - omegaNu - omegaB - omegaC;
  annotation(Documentation(info="<html>
<p><code>&Omega;<sub>&gamma;</sub>h&sup2; = 2.47282&times;10<sup>&minus;5</sup>(T/2.7255 K)<sup>4</sup></code> from the
blackbody energy density; each neutrino species carries <code>(7/8)(4/11)<sup>4/3</sup></code> of the photon
density. <code>H<sub>0</sub>/c</code> is returned in <code>1/Mpc</code>, the unit of every conformal
quantity in this package.</p>
</html>"));
end densities;
