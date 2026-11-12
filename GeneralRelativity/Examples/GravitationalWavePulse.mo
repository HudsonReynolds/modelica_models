within GeneralRelativity.Examples;
model GravitationalWavePulse "A linearized perturbation pulse splits, propagates at c and leaves through the outflow boundaries"
  extends GeneralRelativity.Numerical.LinearizedWave1D(
    n = 201, length = 100.0, c = 1.0, amplitude = 1.0e-3, centre = 50.0, width = 4.0, sourceAmplitude = 0.0);
equation
  assert(hMax <= 1.001 * amplitude + 1.0e-12, "The pulse grew; the discretization is unstable");
  assert(energy <= 1.01 * initialEnergy + 1.0e-18, "Wave energy increased; outflow boundaries are not absorbing");
  assert(time < 80.0 or energy < 0.02 * initialEnergy, "The pulse did not leave the domain by t = 80");
  annotation(experiment(StartTime = 0.0, StopTime = 100.0, Tolerance = 1.0e-8, Interval = 0.5),
    Documentation(info="<html>
<p>The initial Gaussian splits into two halves moving at <code>c</code> in each
direction. By <code>t = 50</code> they reach the boundaries, which absorb them
through the Sommerfeld conditions; by <code>t = 80</code> less than two percent of
the initial energy remains. The amplitude never exceeds its initial value.</p>
</html>"));
end GravitationalWavePulse;
