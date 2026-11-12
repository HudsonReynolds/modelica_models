within GeneralRelativity.Visualization;
impure function writeFieldSlice "Append a scalar field sampled on an axis-aligned plane"
  input String fileName;
  input String name "Field label shown by the viewer";
  input String plane "\"xy\", \"xz\" or \"yz\"";
  input Real offset "Coordinate of the plane along its normal";
  input Real u[:] "First in-plane coordinate samples";
  input Real v[:] "Second in-plane coordinate samples";
  input Real values[size(u, 1), size(v, 1)] "Field values at (u[i], v[j])";
protected
  String text;
algorithm
  text := "{\"type\":\"field\",\"name\":\"" + name + "\",\"plane\":\"" + plane + "\",\"offset\":" + jsonNumber(offset)
    + ",\"u\":" + jsonVector(u) + ",\"v\":" + jsonVector(v) + ",\"values\":[";
  for i in 1:size(u, 1) loop
    text := text + (if i > 1 then "," else "") + jsonVector(values[i, :]);
  end for;
  text := text + "]}";
  Modelica.Utilities.Streams.print(text, fileName);
end writeFieldSlice;
