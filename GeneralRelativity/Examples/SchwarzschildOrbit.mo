within GeneralRelativity.Examples;
model SchwarzschildOrbit "Precessing timelike orbit about a Schwarzschild black hole in geometric units, with curvature checks"
  extends GeneralRelativity.Geodesics.Geodesic(
    redeclare function metric = GeneralRelativity.Metrics.schwarzschild,
    p = {rs},
    c = 1.0,
    x_start = {0.0, r0, GeneralRelativity.Constants.pi / 2.0, 0.0},
    u_start = {E / (1.0 - rs / r0), 0.0, 0.0, L / r0^2},
    step = {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3},
    massive = true);
  parameter Real M = 1.0 "Black hole mass as a length, G = c = 1";
  final parameter Real rs = 2.0 * M "Schwarzschild radius";
  parameter Real r0 = 10.0 * M "Initial radius; a periastron because u^r = 0 and L exceeds the circular value";
  parameter Real L = 4.2 * M "Specific angular momentum; stable orbits need L > 2 sqrt(3) M";
  final parameter Real E = sqrt((1.0 - rs / r0) * (1.0 + L^2 / r0^2)) "Specific energy for a turning point at r0";
  final parameter Real weakFieldAdvanceEstimate = 6.0 * GeneralRelativity.Constants.pi * M / (L^2 / M)
    "Leading-order periastron advance per orbit, 6 pi M / p with p = L^2/M";
  GeneralRelativity.Tensors.Curvature k "Curvature at the particle";
  Real energy "Killing energy -g(d_t, u); constant along the geodesic";
  Real angularMomentum "Killing angular momentum g(d_phi, u); constant along the geodesic";
  Real kretschmannExpected "12 rs^2 / r^6";
  Real ricciNorm "Frobenius norm of the Ricci tensor; zero in vacuum";
  Real properTimeRate "d tau / d t = c / u^t";
  discrete Real periastronPhi(start = 0.0, fixed = true) "Azimuth at the last periastron";
  discrete Real periastronAdvance(start = 0.0, fixed = true) "Azimuth gained per orbit beyond 2 pi, measured between periastra";
  discrete Integer periastra(start = 0, fixed = true) "Periastron passages seen";
equation
  k = GeneralRelativity.Metrics.curvatureAt(metric, x, p, step);
  energy = -GeneralRelativity.Geodesics.killingConserved(g, {1.0, 0.0, 0.0, 0.0}, u);
  angularMomentum = GeneralRelativity.Geodesics.killingConserved(g, {0.0, 0.0, 0.0, 1.0}, u);
  kretschmannExpected = 12.0 * rs^2 / x[2]^6;
  ricciNorm = sqrt(sum(k.ricci[a, b]^2 for a in 1:4, b in 1:4));
  properTimeRate = c / u[1];
  when u[2] > 0.0 then
    periastronPhi = x[4];
    periastra = pre(periastra) + 1;
    periastronAdvance = if pre(periastra) >= 1 then x[4] - pre(periastronPhi) - 2.0 * GeneralRelativity.Constants.pi else 0.0;
  end when;
  assert(x[2] > rs, "The orbit crossed the horizon; the Schwarzschild chart is not valid there");
  assert(abs(energy - E) < 1.0e-5 * E, "Killing energy drifted; tighten the solver tolerance");
  assert(abs(angularMomentum - L) < 1.0e-5 * L, "Killing angular momentum drifted; tighten the solver tolerance");
  assert(abs(normResidual) < 1.0e-5, "Four-velocity normalization drifted");
  assert(abs(k.kretschmann - kretschmannExpected) < 1.0e-4 * kretschmannExpected,
    "Numerical Kretschmann scalar disagrees with 12 rs^2 / r^6");
  assert(ricciNorm < 1.0e-5, "Ricci tensor is not zero in the Schwarzschild vacuum");
  annotation(experiment(StartTime = 0.0, StopTime = 1500.0, Tolerance = 1.0e-10, Interval = 0.5),
    Documentation(info="<html>
<p>A test particle starts at periastron <code>r<sub>0</sub> = 10M</code> with angular
momentum <code>L = 4.2M</code>, which makes the apoapsis about <code>21M</code>. The
affine parameter (Modelica <code>time</code>) is proper time; <code>t</code> is
coordinate time. Three things are checked continuously:</p>
<ul>
<li>the Killing constants <code>E</code> and <code>L</code> and the normalization
<code>g(u,u) = -1</code> stay at their initial values, which tests the geodesic
integration and the finite-difference Christoffel symbols;</li>
<li>the Ricci tensor computed at the particle from second differences of the
metric vanishes and the Kretschmann scalar equals
<code>12r<sub>s</sub>&sup2;/r<sup>6</sup></code>, which tests the curvature pipeline;</li>
<li>the chart stays outside the horizon.</li>
</ul>
<p><code>periastronAdvance</code> is the azimuth gained beyond <code>2&pi;</code>
between successive periastra: about 2.14 rad per orbit here, twice the
leading-order <code>6&pi;M/p = 1.07</code> rad estimate, because the orbit is deep in
the strong field; the two agree only when <code>r<sub>0</sub> &gt;&gt; r<sub>s</sub></code>.</p>
<p>This model has no file output, so it builds on every OpenModelica target. <a href=\"modelica://GeneralRelativity.Examples.SchwarzschildOrbitScene\">SchwarzschildOrbitScene</a> adds the Three.js scene stream, which needs <code>Modelica.Utilities</code> to link.</p>
</html>"));
end SchwarzschildOrbit;
