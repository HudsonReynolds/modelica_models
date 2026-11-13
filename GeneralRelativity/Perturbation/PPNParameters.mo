within GeneralRelativity.Perturbation;
record PPNParameters "The ten parametrized post-Newtonian parameters; defaults are the general-relativity values"
  Real gamma = 1.0 "Space curvature produced by unit rest mass";
  Real beta = 1.0 "Nonlinearity in the superposition law for gravity";
  Real xi = 0.0 "Preferred-location effects";
  Real alpha1 = 0.0 "Preferred-frame effect";
  Real alpha2 = 0.0 "Preferred-frame effect";
  Real alpha3 = 0.0 "Preferred-frame effect and violation of momentum conservation";
  Real zeta1 = 0.0 "Violation of total momentum conservation";
  Real zeta2 = 0.0 "Violation of total momentum conservation";
  Real zeta3 = 0.0 "Violation of total momentum conservation";
  Real zeta4 = 0.0 "Violation of total momentum conservation";
  annotation(Documentation(info="<html>
<p>Will's PPN parameters (<i>Theory and Experiment in Gravitational Physics</i>,
Living Reviews in Relativity 17:4). General relativity has
<code>&gamma; = &beta; = 1</code> and all others zero. Fully conservative theories have
<code>&alpha;<sub>i</sub> = &zeta;<sub>i</sub> = 0</code>; semi-conservative ones have
<code>&zeta;<sub>i</sub> = &alpha;<sub>3</sub> = 0</code>. Observational bounds (Cassini,
lunar laser ranging, pulsars) are <code>|&gamma; - 1| &lt; 2.3&times;10<sup>-5</sup></code> and
<code>|&beta; - 1| &lt; 8&times;10<sup>-5</sup></code>. See
<a href=\"modelica://GeneralRelativity.Perturbation.bransDicke\">bransDicke</a> and
<a href=\"modelica://GeneralRelativity.Perturbation.nordtvedtEta\">nordtvedtEta</a>.</p>
</html>"));
end PPNParameters;
