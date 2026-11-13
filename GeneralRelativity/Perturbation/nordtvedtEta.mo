within GeneralRelativity.Perturbation;
function nordtvedtEta "Nordtvedt parameter eta_N of the strong-equivalence-principle violation"
  input PPNParameters parameters "PPN parameters";
  output Real eta "eta_N = 4 beta - gamma - 3 - 10/3 xi - alpha1 + 2/3 alpha2 - 2/3 zeta1 - 1/3 zeta2";
algorithm
  eta := 4.0 * parameters.beta - parameters.gamma - 3.0 - 10.0 / 3.0 * parameters.xi - parameters.alpha1
    + 2.0 / 3.0 * parameters.alpha2 - 2.0 / 3.0 * parameters.zeta1 - 1.0 / 3.0 * parameters.zeta2;
  annotation(Documentation(info="<html>
<p>A body with gravitational self-energy <code>&Omega;</code> and rest energy
<code>Mc&sup2;</code> has passive-to-inertial mass ratio
<code>1 + &eta;<sub>N</sub>&Omega;/(Mc&sup2;)</code> (Will 2014). It vanishes in general
relativity; lunar laser ranging bounds it to about <code>4&times;10<sup>-4</sup></code>.</p>
</html>"));
end nordtvedtEta;
