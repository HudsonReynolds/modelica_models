within GeneralRelativity.Visualization;
record SpacetimePoint "Observer-level quantities at one event of a 2-D slice, the data behind the river, clock and light-cone styles"
  Real valid "1 where Eulerian observers exist (g^00 < 0), else 0; a Real so the record can sit in equations";
  Real lapse "alpha = 1/sqrt(-g^00): proper-time rate of Eulerian (zero-angular-momentum) observers per unit coordinate time";
  Real timeLag "1 - alpha, computed without cancellation so weak fields keep their digits";
  Real riverSpeed "Speed of the river of space, sqrt(1 + g_00), in units of c: free fall from rest at infinity seen by static observers; exceeds 1 inside the static limit";
  Real dragSpeed "Frame-dragging flow -g(e_t, e_phi)/g(e_phi, e_phi) along the display tangential direction, display length per unit ct";
  Real coneValid "1 where both light-cone sections below exist";
  Real coneOutward "Outward radial coordinate speed of light, display length per unit ct";
  Real coneInward "Inward radial coordinate speed of light (negative outside a horizon)";
  Real conePrograde "Tangential coordinate speed of light along +phi";
  Real coneRetrograde "Tangential coordinate speed of light along -phi (negative)";
  annotation(Documentation(info="<html>
<p>Filled by <a href=\"modelica://GeneralRelativity.Visualization.spacetimePointAt\">spacetimePointAt</a>.
Speeds are measured in the <i>display</i> directions of the slice: radial, and tangential
arc length <code>R d&phi;</code> in the plane, per unit <code>x<sup>0</sup> = ct</code>, so a light cone far
from any mass has speeds <code>&plusmn;1</code>.</p>
</html>"));
end SpacetimePoint;
