within GeneralRelativity.Examples;
model IndexNotationDemo "Tensor expressions typed as text on the Schwarzschild curvature: metric trace, Kretschmann scalar, Ricci contraction, Einstein identity, Levi-Civita square"
  parameter Real rs = 2.0 "Schwarzschild radius in geometric units";
  parameter Real r = 10.0 "Sample radius";
  parameter Real step[4] = {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3} "Finite-difference steps";
  final parameter Real values[6] = report(rs, r, step);
  Real metricTrace = values[1] "g_{mu nu} g^{mu nu}, equal to 4";
  Real kretschmann = values[2] "R_{abcd} R^{abcd}";
  Real kretschmannClosedForm = 12.0 * rs^2 / r^6 "12 rs^2 / r^6";
  Real ricciResidual = values[3] "Largest |R^a_{bac}|, zero in vacuum";
  Real einsteinResidual = values[4] "Largest |R_{ab} - 1/2 g_{ab} R - G_{ab}|";
  Real leviCivitaSquare = values[5] "epsilon_{abcd} epsilon^{abcd}, equal to -24";
  Real staticEnergyDensity = values[6] "T_{ab} u^a u^b for the static observer, zero in vacuum";
initial algorithm
  printExpressions(rs, r, step);
equation
  assert(abs(metricTrace - 4.0) < 1.0e-12, "g_{mu nu} g^{mu nu} is not 4");
  assert(abs(kretschmann / kretschmannClosedForm - 1.0) < 1.0e-4, "R_{abcd} R^{abcd} disagrees with 12 rs^2 / r^6");
  assert(ricciResidual < 1.0e-5, "R^a_{bac} is not zero in vacuum");
  assert(einsteinResidual < 1.0e-12, "R_{ab} - 1/2 g_{ab} R - G_{ab} is not zero");
  assert(abs(leviCivitaSquare + 24.0) < 1.0e-10, "epsilon_{abcd} epsilon^{abcd} is not -24");
  assert(abs(staticEnergyDensity) < 1.0e-5, "T_{ab} u^a u^b is not zero in vacuum");
protected
  function system "Schwarzschild curvature at the sample radius with the static observer's four-velocity u"
    input Real rs;
    input Real r;
    input Real step[4];
    output GeneralRelativity.IndexNotation.TensorSystem sys;
  algorithm
    sys := GeneralRelativity.IndexNotation.systemAt(GeneralRelativity.Metrics.schwarzschild,
      {0.0, r, GeneralRelativity.Constants.pi / 2.0, 0.0}, {rs}, step);
    sys := GeneralRelativity.IndexNotation.defineVector(sys, "u", {1.0 / sqrt(1.0 - rs / r), 0.0, 0.0, 0.0});
  end system;

  function report "Evaluate the expressions and return the numbers"
    input Real rs;
    input Real r;
    input Real step[4];
    output Real values[6];
  protected
    import GeneralRelativity.IndexNotation.scalar;
    import GeneralRelativity.IndexNotation.matrix;
    GeneralRelativity.IndexNotation.TensorSystem sys = system(rs, r, step);
    Real M[4, 4];
  algorithm
    values[1] := scalar(sys, "g_{mu nu} g^{mu nu}");
    values[2] := scalar(sys, "R_{abcd} R^{abcd}");
    M := matrix(sys, "R^{a}_{bac}");
    values[3] := max(abs(M[i, j]) for i in 1:4, j in 1:4);
    M := matrix(sys, "R_{ab} - 1/2 g_{ab} R - G_{ab}");
    values[4] := max(abs(M[i, j]) for i in 1:4, j in 1:4);
    values[5] := scalar(sys, "epsilon_{abcd} epsilon^{abcd}");
    values[6] := scalar(sys, "T_{ab} u^a u^b");
  end report;

  impure function printExpressions "Print a few expressions and their values; edit the strings to try your own"
    input Real rs;
    input Real r;
    input Real step[4];
  protected
    import GeneralRelativity.IndexNotation.evaluate;
    import GeneralRelativity.IndexNotation.toString;
    import Modelica.Utilities.Streams.print;
    GeneralRelativity.IndexNotation.TensorSystem sys = system(rs, r, step);
    String expressions[8] = {
      "g_{mu nu} g^{mu nu}",
      "R_{abcd} R^{abcd}",
      "K",
      "R^{a}_{bac}",
      "R_{ab} - 1/2 g_{ab} R - G_{ab}",
      "epsilon_{abcd} epsilon^{abcd}",
      "Gamma^i_{jk}",
      "g_{ab} u^a u^b"};
  algorithm
    for i in 1:size(expressions, 1) loop
      print(expressions[i] + "  ->  " + toString(evaluate(sys, expressions[i])));
    end for;
  end printExpressions;
  annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 1.0e-8, Interval = 0.1),
    Documentation(info="<html>
<p>Every quantity here is typed as a string and evaluated by
<a href=\"modelica://GeneralRelativity.IndexNotation\">IndexNotation</a> against
the Schwarzschild curvature at <code>r = 5 r<sub>s</sub></code>. The expressions and
their values are printed at initialization, which is the quickest way to try
your own: edit the strings in <code>printExpressions</code> and simulate.</p>
</html>"));
end IndexNotationDemo;
