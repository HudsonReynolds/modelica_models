within GeneralRelativity.Examples;
model LightDeflection "A light ray passing a Schwarzschild mass: the null geodesic's bending angle against 4GM/(c^2 b) and its higher-order corrections"
  extends GeneralRelativity.Geodesics.Geodesic(
    redeclare function metric = GeneralRelativity.Metrics.schwarzschild,
    p = {2.0 * M},
    c = 1.0,
    x_start = {0.0, r0, GeneralRelativity.Constants.pi / 2.0, phi0},
    u_start = {1.0 / (1.0 - 2.0 * M / r0), -sqrt(1.0 - (1.0 - 2.0 * M / r0) * b^2 / r0^2), 0.0, b / r0^2},
    step = {1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4},
    massive = false);
  // Lengths in units of the impact parameter: the solver's absolute tolerance then resolves u^phi ~ b/r^2
  // far out, which is tiny when lengths are in units of M.
  parameter Real b = 1.0 "Impact parameter L/E (the length unit)";
  parameter Real M = b / 200.0 "Mass as a length, G = c = 1";
  parameter Real r0 = 200.0 * b "Start and end radius; the azimuth swept beyond it is added analytically";
  final parameter Real phi0 = GeneralRelativity.Constants.pi - asin(b / r0) "Start on the far left of the hole, moving toward +x";
  final parameter Real m = M / b "Expansion parameter";
  final parameter Real deflectionSeries = 4.0 * m + 15.0 * GeneralRelativity.Constants.pi / 4.0 * m^2 + 128.0 / 3.0 * m^3
    + 3465.0 * GeneralRelativity.Constants.pi / 64.0 * m^4 "Bending angle to fourth order in M/b (Keeton and Petters 2005)";
  final parameter Real deflectionWeakField = 4.0 * m "Thin-lens (weak-field) bending angle 4GM/(c^2 b) used by Lensing";
  Real rayEnergy "Killing energy -g(d_t, u); 1 here";
  Real rayAngularMomentum "g(d_phi, u) = b";
  discrete Real deflection(start = 0.0, fixed = true) "Measured bending angle";
equation
  rayEnergy = -GeneralRelativity.Geodesics.killingConserved(g, {1.0, 0.0, 0.0, 0.0}, u);
  rayAngularMomentum = GeneralRelativity.Geodesics.killingConserved(g, {0.0, 0.0, 0.0, 1.0}, u);
  assert(abs(normResidual) < 1.0e-8, "The ray is no longer null");
  // A sanity check only: far out u^phi ~ b/r^2 is near the solver's absolute tolerance, where it barely
  // contributes to the bending; the deflection itself is asserted to 1e-6 below.
  assert(abs(rayAngularMomentum / rayEnergy - b) < 1.0e-5 * b, "The impact parameter L/E drifted");
  when u[2] > 0.0 and x[2] >= r0 then
    // Bending = total swept azimuth - pi. Beyond r0 on each side the ray sweeps
    // int_r0^inf (b/r^2)/sqrt(1 - (1 - 2M/r) b^2/r^2) dr = asin(b/r0) - M b^3/(4 r0^4) + ..., so adding 2 asin(b/r0)
    // is exact to 1e-11 here.
    deflection = x[4] - phi0 + 2.0 * asin(b / r0) - GeneralRelativity.Constants.pi;
    assert(abs(deflection - deflectionSeries) < 1.0e-6 * deflectionSeries,
      "Null-geodesic deflection disagrees with the fourth-order series");
    assert(abs(deflection / deflectionWeakField - 1.0) < 2.0 * m * 15.0 * GeneralRelativity.Constants.pi / 16.0,
      "Deflection is not 4GM/(c^2 b) to first order");
    terminate("The ray is back at r0");
  end when;
  annotation(experiment(StartTime = 0.0, StopTime = 2000.0, Tolerance = 1.0e-12, Interval = 0.5),
    Documentation(info="<html>
<p>A light ray with impact parameter <code>b = 200M</code> is integrated as a null geodesic of the
Schwarzschild metric from <code>r = 200b</code>, past the mass and back out to the same radius
(lengths are in units of <code>b</code>, which keeps the integration well conditioned). The
bending angle (swept azimuth minus &pi;, with the azimuth swept beyond <code>r<sub>0</sub></code>,
<code>asin(b/r<sub>0</sub>)</code> on each side up to <code>O(Mb&sup3;/r<sub>0</sub><sup>4</sup>)</code>, added back) is
asserted against the series <code>4m + (15&pi;/4)m&sup2; + (128/3)m&sup3; + (3465&pi;/64)m<sup>4</sup></code>,
<code>m = M/b</code>, to one part in 10<sup>6</sup>. Its leading term is the thin-lens deflection
<code>4GM/(c&sup2;b)</code> that <a href=\"modelica://GeneralRelativity.Lensing\">Lensing</a> uses, so this model ties
the lensing package to the geodesic equation.</p>
</html>"));
end LightDeflection;
