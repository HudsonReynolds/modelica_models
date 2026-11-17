within GeneralRelativity.CompactBinaries;
function circularOrbitSpeed "Relative speed of an exactly circular orbit of the conservative post-Newtonian dynamics at separation r"
  input Real r "Harmonic-coordinate separation";
  input Real M = 1.0 "Total mass";
  input Real eta = 0.25 "Symmetric mass ratio";
  input Integer pnOrder = 2 "0, 1 or 2";
  output Real speed "|v| with v perpendicular to x";
protected
  Real a[3];
algorithm
  // Fixed point of v^2 = r |a_radial(v)| for the conservative acceleration (no radiation reaction).
  speed := sqrt(M / r);
  for i in 1:60 loop
    a := relativeAcceleration({r, 0.0, 0.0}, {0.0, speed, 0.0}, M, eta, pnOrder, false);
    speed := sqrt(max(-a[1] * r, 0.0));
  end for;
  annotation(Documentation(info="<html>
<p>The speed for which the centripetal acceleration <code>v&sup2;/r</code> equals the radial pull of
<a href=\"modelica://GeneralRelativity.CompactBinaries.relativeAcceleration\">relativeAcceleration</a> (conservative
terms, <code>&#7563; = 0</code>), found by fixed-point iteration. It agrees with
<code>r</code>&middot;<a href=\"modelica://GeneralRelativity.CompactBinaries.circularOmega\">circularOmega</a> up to the
truncation of that series (relative <code>O(&gamma;&sup3;)</code>), but is exactly circular for the equations of motion:
starting a 2PN orbit with the series value leaves an eccentricity of about 0.01 at <code>r = 14M</code>.</p>
</html>"));
end circularOrbitSpeed;
