within GeneralRelativity.Visualization;
impure function writeView "Choose the viewer style a scene opens with and preset its options"
  input String fileName;
  input String style = "embedding" "embedding, river, clocks, spacetime or raytrace";
  input String optionNames[:] = fill("", 0) "Style option identifiers, as listed in the viewer README";
  input Real optionValues[size(optionNames, 1)] = fill(0.0, size(optionNames, 1)) "Option values; booleans as 0 or 1";
protected
  String options;
algorithm
  options := "{";
  for i in 1:size(optionNames, 1) loop
    options := options + (if i > 1 then "," else "") + "\"" + optionNames[i] + "\":" + jsonNumber(optionValues[i]);
  end for;
  Modelica.Utilities.Streams.print("{\"type\":\"view\",\"style\":\"" + style + "\",\"options\":" + options + "}}", fileName);
  annotation(Documentation(info="<html>
<p>Only a starting point: every style whose data is in the scene can be selected in the viewer
at run time, and its options changed there.</p>
</html>"));
end writeView;
