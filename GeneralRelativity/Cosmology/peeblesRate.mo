within GeneralRelativity.Cosmology;
function peeblesRate "dx/dz of the hydrogen ionization fraction in Peebles' three-level atom (RECFAST fudge 1.14)"
  input Real z "Redshift";
  input Real x "Ionized fraction of hydrogen, n_p / n_H";
  input Real xHe "Electrons per hydrogen atom from helium";
  input Parameters par;
  output Real dxdz;
protected
  constant Real kB = 1.380649e-23;
  constant Real hP = 6.62607015e-34;
  constant Real me = 9.1093837015e-31;
  constant Real eV = 1.602176634e-19;
  constant Real B1 = 13.605693 * eV "Hydrogen ground-state binding energy";
  constant Real B2 = B1 / 4.0 "Binding energy of n = 2";
  constant Real lambdaAlpha = 121.5670e-9 "Lyman-alpha wavelength";
  constant Real Lambda2s = 8.22458 "Two-photon decay rate of 2s, 1/s";
  Real T = par.TCMB * (1.0 + z) "Matter temperature = radiation temperature";
  Real nH = hydrogenDensity(par) * (1.0 + z)^3;
  Real H = hubble(1.0 / (1.0 + z), par) * 1000.0 / GeneralRelativity.Constants.Mpc "1/s";
  Real t4 = T / 1.0e4;
  Real alphaB = 1.14 * 1.0e-19 * 4.309 * t4^(-0.6166) / (1.0 + 0.6703 * t4^0.5300) "Case-B recombination, m^3/s";
  Real betaB = alphaB * (2.0 * GeneralRelativity.Constants.pi * me * kB * T / hP^2)^1.5 * exp(-B2 / (kB * T)) "Photoionization from n = 2, 1/s";
  Real n1s = max((1.0 - x) * nH, 1.0e-30);
  Real LambdaAlpha = 8.0 * GeneralRelativity.Constants.pi * H / (lambdaAlpha^3 * n1s) "Escape rate of Lyman-alpha photons";
  Real C = (Lambda2s + LambdaAlpha) / (Lambda2s + LambdaAlpha + betaB) "Peebles factor: probability that n = 2 reaches the ground state";
  Real dxdt;
algorithm
  dxdt := -C * (alphaB * (x + xHe) * nH * x - betaB * (1.0 - x) * exp(-(B1 - B2) / (kB * T)));
  dxdz := -dxdt / ((1.0 + z) * H);
  annotation(Documentation(info="<html>
<p>P. J. E. Peebles, ApJ 153, 1 (1968): recombinations go through the <code>n = 2</code> level, which
reaches the ground state only by the two-photon decay of <code>2s</code> or by Lyman-&alpha; photons
redshifting out of the line, so</p>
<p><code>dx/dt = &minus;C [&alpha;<sub>B</sub> n<sub>e</sub> x &minus; &beta;<sub>B</sub>(1 &minus; x) e<sup>&minus;h&nu;<sub>&alpha;</sub>/kT</sup>]</code>,
<code>C = (&Lambda;<sub>2s</sub> + &Lambda;<sub>&alpha;</sub>)/(&Lambda;<sub>2s</sub> + &Lambda;<sub>&alpha;</sub> + &beta;<sub>B</sub>)</code>,
<code>&Lambda;<sub>&alpha;</sub> = 8&pi;H/(&lambda;<sub>&alpha;</sub>&sup3;n<sub>1s</sub>)</code>.</p>
<p>The case-B coefficient is the P&eacute;quignot et al. (1991) fit times the RECFAST fudge factor 1.14
(Seager, Sasselov &amp; Scott 1999), which reproduces the multilevel atom to about 1%; detailed balance gives
<code>&beta;<sub>B</sub></code>. The matter temperature is taken equal to the radiation temperature (true to
<code>z &asymp; 200</code>).</p>
</html>"));
end peeblesRate;
