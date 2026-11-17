within Tests;
model GeneralRelativityBinaryInspiralTest "A Newtonian binary with 2.5PN radiation reaction shrinks on Peters' schedule"
  parameter Real eta = 0.25 "Symmetric mass ratio (equal masses)";
  parameter Real r0 = 10.0 "Initial separation in M";
  parameter Real rEnd = 4.0 "Separation at which the time is checked";
  final parameter Real tExpected = 5.0 / 256.0 * (r0^4 - rEnd^4) / eta "Peters: (5/256)(r0^4 - r^4)/(eta M^3)";
  Real x[3](start = {r0, 0.0, 0.0}, each fixed = true);
  Real v[3](start = {GeneralRelativity.CompactBinaries.inspiralRadialVelocity(r0, 1.0, eta), sqrt(1.0 / r0), 0.0}, each fixed = true);
  Real separation;
equation
  der(x) = v;
  der(v) = GeneralRelativity.CompactBinaries.relativeAcceleration(x, v, 1.0, eta, 0, true);
  separation = sqrt(x * x);
  assert(time < 1.5 * tExpected, "The binary has not shrunk to rEnd in time");
  when separation <= rEnd then
    assert(abs(time / tExpected - 1.0) < 0.01, "Inspiral time differs from Peters' formula by more than 1%");
    terminate("Reached rEnd on Peters' schedule");
  end when;
  annotation(experiment(StartTime = 0.0, StopTime = 2000.0, Tolerance = 1.0e-10, Interval = 0.5),
    Documentation(info="<html>
<p>Newtonian gravity plus the 2.5PN radiation reaction of
<a href=\"modelica://GeneralRelativity.CompactBinaries.relativeAcceleration\">CompactBinaries.relativeAcceleration</a>,
started quasi-circular (Kepler speed and the leading-order inward drift) at <code>r<sub>0</sub> = 10M</code>,
must reach <code>4M</code> after <code>(5/256)(r<sub>0</sub>&#8308; &minus; r&#8308;)/(&eta;M&sup3;)</code> = 761 M within 1%.</p>
</html>"));
end GeneralRelativityBinaryInspiralTest;
