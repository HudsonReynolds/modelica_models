within GeneralRelativity.Visualization;
function spacetimePointAt "Lapse, river speed, frame dragging and light-cone speeds at an event in the equatorial plane of any metric"
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric function";
  input Real x[4] "Event: {ct, r, pi/2, phi} in a spherical chart or {ct, x, y, 0} in a Cartesian chart";
  input Real p[:] "Metric parameters";
  input Boolean spherical = true "true for {ct, r, theta, phi}, false for {ct, x, y, z}";
  output SpacetimePoint s;
protected
  Real g[4, 4];
  Real gInv[4, 4];
  Boolean accepted;
  Real eTime[4] = {1.0, 0.0, 0.0, 0.0};
  Real eRadial[4];
  Real eTangential[4];
  Real planeRadius;
  Real q;
  Real root;
  Real a;
  Real b;
  Real discriminant;
  Boolean radialOk;
  Boolean tangentialOk;
algorithm
  g := metric(x, p);
  (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(g);
  // Display directions: unit radial and unit arc length in the plane.
  if spherical then
    planeRadius := x[2] * sin(x[3]);
    eRadial := {0.0, 1.0, 0.0, 0.0};
    eTangential := {0.0, 0.0, 0.0, 1.0 / planeRadius};
  else
    planeRadius := sqrt(x[2]^2 + x[3]^2);
    eRadial := {0.0, x[2] / planeRadius, x[3] / planeRadius, 0.0};
    eTangential := {0.0, -x[3] / planeRadius, x[2] / planeRadius, 0.0};
  end if;
  q := -gInv[1, 1];
  s.valid := if accepted and q > 0.0 then 1.0 else 0.0;
  if s.valid > 0.5 then
    root := sqrt(q);
    s.lapse := 1.0 / root;
    // 1 - 1/sqrt(q) = (q - 1) / (sqrt(q) (sqrt(q) + 1)): no cancellation when q is close to 1.
    s.timeLag := (q - 1.0) / (root * (root + 1.0));
  else
    s.lapse := 0.0;
    s.timeLag := 1.0;
  end if;
  s.riverSpeed := sqrt(max(1.0 + g[1, 1], 0.0));
  s.dragSpeed := -(eTime * (g * eTangential)) / (eTangential * (g * eTangential));
  // Light cone in the (t, e) plane: g(e_t + v e, e_t + v e) = 0, i.e. a v^2 + 2 b v + g_tt = 0.
  a := eRadial * (g * eRadial);
  b := eTime * (g * eRadial);
  discriminant := b^2 - a * g[1, 1];
  radialOk := a > 0.0 and discriminant >= 0.0;
  s.coneOutward := if radialOk then (-b + sqrt(discriminant)) / a else 0.0;
  s.coneInward := if radialOk then (-b - sqrt(discriminant)) / a else 0.0;
  a := eTangential * (g * eTangential);
  b := eTime * (g * eTangential);
  discriminant := b^2 - a * g[1, 1];
  tangentialOk := a > 0.0 and discriminant >= 0.0;
  s.conePrograde := if tangentialOk then (-b + sqrt(discriminant)) / a else 0.0;
  s.coneRetrograde := if tangentialOk then (-b - sqrt(discriminant)) / a else 0.0;
  s.coneValid := if radialOk and tangentialOk then 1.0 else 0.0;
  annotation(Documentation(info="<html>
<p>All quantities are evaluated from the metric alone, so every
<a href=\"modelica://GeneralRelativity.Interfaces.MetricFunction\">MetricFunction</a> gets the
river, clock and light-cone visualizations without extra code.</p>
<ul>
<li><b>Lapse</b> <code>&alpha; = 1/&radic;(&minus;g<sup>00</sup>)</code> is the rate of the clocks of observers at
rest in the slicing (static observers for Schwarzschild; zero-angular-momentum observers for Kerr).</li>
<li><b>River speed</b> <code>v = &radic;(1 + g<sub>00</sub>)</code>: a body dropped from rest at infinity has
Killing energy 1, so a static observer measures Lorentz factor <code>1/&radic;(&minus;g<sub>00</sub>)</code>.
In the river model of Hamilton and Lisle (<i>Am. J. Phys.</i> 76, 519 (2008)) space itself flows inward at
this speed; for Schwarzschild <code>v = &radic;(r<sub>s</sub>/r)</code>, the Newtonian escape speed, reaching
<code>c</code> at the horizon. The formula needs <code>x<sup>0</sup></code> to be the time of the stationary
Killing vector normalized at infinity; for Kerr it reaches <code>c</code> at the static limit, and the
rotation of the flow is <b>frame dragging</b> <code>&minus;g(e<sub>t</sub>, e<sub>&phi;</sub>)/g(e<sub>&phi;</sub>, e<sub>&phi;</sub>)</code>.</li>
<li><b>Light cones</b>: coordinate speeds of light along the display radial and tangential directions,
from <code>g(e<sub>t</sub> + v e, e<sub>t</sub> + v e) = 0</code>. For Schwarzschild coordinates the radial speeds
are <code>&plusmn;(1 &minus; r<sub>s</sub>/r)</code>, so cones close up toward the horizon.</li>
</ul>
<p>In a Cartesian chart the radial direction points away from the chart origin, which should be the
centre of the dominant mass.</p>
</html>"));
end spacetimePointAt;
