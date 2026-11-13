within GeneralRelativity.NewtonCartan;
function fieldEquationResidual "Newton-Cartan field equation R_ab - (4 pi G / c^2) rho tau_a tau_b, which is Poisson's equation"
  input Real Ric[4, 4] "Ricci tensor of the connection, R_ab = R^r_arb";
  input Real rho(unit = "kg/m3") "Mass density";
  input Real gravitationalConstant(unit = "m3/(kg.s2)") = 1.0 "G; 1 in geometric units";
  input Real c(unit = "m/s") = 1.0 "Speed of light; 1 in geometric units";
  output Real residual[4, 4] "Zero where the connection solves the field equation";
protected
  Real tau[4] = clockForm();
algorithm
  for a in 1:4 loop
    for b in 1:4 loop
      residual[a, b] := Ric[a, b] - 4.0 * GeneralRelativity.Constants.pi * gravitationalConstant * rho / c^2 * tau[a] * tau[b];
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>With <code>x<sup>0</sup> = ct</code>, <code>R<sub>00</sub> = &part;<sub>i</sub>&Gamma;<sup>i</sup><sub>00</sub> = &nabla;&sup2;&Phi;/c&sup2;</code>,
so the <code>00</code> component is Poisson's equation <code>&nabla;&sup2;&Phi; = 4&pi;G&rho;</code> and the
others vanish. Compare the general-relativistic
<a href=\"modelica://GeneralRelativity.StressEnergy.einsteinResidual\">StressEnergy.einsteinResidual</a>.</p>
</html>"));
end fieldEquationResidual;
