within GeneralRelativity.CompactBinaries;
function relativeAcceleration "Post-Newtonian relative acceleration of a compact binary in harmonic coordinates, through 2PN (conservative) and 2.5PN (radiation reaction)"
  input Real x[3] "Relative position x1 - x2";
  input Real v[3] "Relative velocity dx/dt";
  input Real M = 1.0 "Total mass m1 + m2 (G = c = 1)";
  input Real eta = 0.25 "Symmetric mass ratio m1 m2 / M^2";
  input Integer pnOrder = 2 "Conservative terms: 0 Newtonian, 1 first post-Newtonian, 2 second post-Newtonian";
  input Boolean radiationReaction = true "Add the 2.5PN radiation-reaction terms";
  output Real a[3] "Relative acceleration d2x/dt2";
protected
  Real r = sqrt(x * x);
  Real n[3] = x / r;
  Real rdot = n * v;
  Real v2 = v * v;
  Real gamma = M / r;
  Real A = 0.0;
  Real B = 0.0;
algorithm
  if pnOrder >= 1 then
    A := A - 1.5 * rdot^2 * eta + v2 + 3.0 * eta * v2 - gamma * (4.0 + 2.0 * eta);
    B := B - 4.0 * rdot + 2.0 * rdot * eta;
  end if;
  if pnOrder >= 2 then
    A := A + 15.0 / 8.0 * rdot^4 * eta - 45.0 / 8.0 * rdot^4 * eta^2 - 4.5 * rdot^2 * eta * v2 + 6.0 * rdot^2 * eta^2 * v2
      + 3.0 * eta * v2^2 - 4.0 * eta^2 * v2^2
      + gamma * (-2.0 * rdot^2 - 25.0 * rdot^2 * eta - 2.0 * rdot^2 * eta^2 - 6.5 * eta * v2 + 2.0 * eta^2 * v2)
      + gamma^2 * (9.0 + 87.0 / 4.0 * eta);
    B := B + 4.5 * rdot^3 * eta + 3.0 * rdot^3 * eta^2 - 7.5 * rdot * eta * v2 - 2.0 * rdot * eta^2 * v2
      + gamma * (2.0 * rdot + 20.5 * rdot * eta + 4.0 * rdot * eta^2);
  end if;
  if radiationReaction then
    A := A - 24.0 / 5.0 * rdot * eta * v2 * gamma - 136.0 / 15.0 * rdot * eta * gamma^2;
    B := B + 8.0 / 5.0 * eta * v2 * gamma + 24.0 / 5.0 * eta * gamma^2;
  end if;
  a := -(M / r^2) * ((1.0 + A) * n + B * v);
  annotation(Documentation(info="<html>
<p><code>a = &minus;(M/r&sup2;)[(1 + A) n + B v]</code> with <code>n = x/r</code>, <code>&#7563; = n&middot;v</code> and
<code>&gamma; = M/r</code>; the coefficients <code>A</code>, <code>B</code> are the harmonic-gauge expressions of
L. Blanchet, <i>Living Rev. Relativ.</i> 17, 2 (2014), section 7 (relative acceleration of the
center-of-mass frame): 1PN and 2PN conservative terms, and the 2.5PN terms
<code>A<sub>2.5</sub> = &minus;(24/5)&eta;&#7563;v&sup2;&gamma; &minus; (136/15)&eta;&#7563;&gamma;&sup2;</code>,
<code>B<sub>2.5</sub> = (8/5)&eta;v&sup2;&gamma; + (24/5)&eta;&gamma;&sup2;</code> of radiation reaction.</p>
<p>Checks in <code>Tests.GeneralRelativityCompactBinaryTests</code>: on a circular orbit the conservative part
gives <a href=\"modelica://GeneralRelativity.CompactBinaries.circularOmega\">circularOmega</a> up to
<code>O(&gamma;&sup3;)</code>, and the radiation-reaction part removes energy at the quadrupole rate
<code>dE/dt = &minus;(32/5)&eta;&sup2;M<sup>5</sup>/r<sup>5</sup></code>.</p>
</html>"));
end relativeAcceleration;
