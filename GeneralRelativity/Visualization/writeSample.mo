within GeneralRelativity.Visualization;
impure function writeSample "Append one trajectory sample: affine parameter, coordinates and attached scalars"
  input String fileName;
  input Real lambda "Affine parameter (proper time for massive particles)";
  input Real x[4] "Coordinates in the scene chart; x[1] = c*t";
  input Real scalars[:] = fill(0.0, 0) "Values such as the Kretschmann scalar at the particle";
algorithm
  Modelica.Utilities.Streams.print("{\"type\":\"sample\",\"lambda\":" + jsonNumber(lambda) + ",\"x\":" + jsonVector(x)
    + ",\"scalars\":" + jsonVector(scalars) + "}", fileName);
end writeSample;
