within GeneralRelativity.Visualization;
impure function writeSample "Append one trajectory sample: affine parameter, coordinates, attached scalars and the track it belongs to"
  input String fileName;
  input Real lambda "Affine parameter (proper time for massive particles)";
  input Real x[4] "Coordinates in the scene chart; x[1] = c*t";
  input Real scalars[:] = fill(0.0, 0) "Values such as the Kretschmann scalar at the particle";
  input String track = "" "Track (body) name declared with writeTrack; empty for the single default track";
algorithm
  Modelica.Utilities.Streams.print("{\"type\":\"sample\",\"lambda\":" + jsonNumber(lambda) + ",\"x\":" + jsonVector(x)
    + ",\"scalars\":" + jsonVector(scalars) + (if track == "" then "" else ",\"track\":\"" + track + "\"") + "}", fileName);
  annotation(Documentation(info="<html>
<p>The viewer animates every track against coordinate time <code>t = x[1]/c</code>, so
bodies integrated with different affine parameters (proper time, or
<code>ct</code> for a Newtonian body) move in step when their samples carry the
same <code>x[1]</code>.</p>
</html>"));
end writeSample;
