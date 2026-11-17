within GeneralRelativity.Visualization;
impure function writeLensImage "Append an image on an angular grid for the viewer's lensfit style: observed, initial or fitted simulation, or residual"
  input String fileName;
  input String name "\"observed\", \"initial\", \"model\" or \"residual\"; other names are shown as extra panels";
  input Real u[:] "First image coordinate of the pixel centres, in arcseconds";
  input Real v[:] "Second image coordinate, in arcseconds";
  input Real values[size(u, 1), size(v, 1)] "Surface brightness at (u[i], v[j])";
  input String unit = "" "Label for the values, for example \"sigma\" for a normalized residual";
protected
  String text;
algorithm
  text := "{\"type\":\"lensImage\",\"name\":\"" + name + "\",\"unit\":\"" + unit + "\",\"u\":" + jsonVector(u, 9) + ",\"v\":" + jsonVector(v, 9) + ",\"values\":[";
  for i in 1:size(u, 1) loop
    text := text + (if i > 1 then "," else "") + jsonVector(values[i, :], 7);
  end for;
  text := text + "]}";
  Modelica.Utilities.Streams.print(text, fileName);
end writeLensImage;
