within GeneralRelativity.Visualization;
impure function writeSliceField "Append one scalar field of a named slice: values[i][j] at (u[i], v[j])"
  input String fileName;
  input String slice "Slice name given to writeSpacetimeSlice";
  input String name "Field name; the viewer styles read lapse, timeLag, river, drag and the cone fields";
  input Real values[:, :] "Field values on the slice grid";
  input Integer significantDigits = 9;
protected
  String text;
algorithm
  text := "{\"type\":\"sliceField\",\"slice\":\"" + slice + "\",\"name\":\"" + name + "\",\"values\":[";
  for i in 1:size(values, 1) loop
    text := text + (if i > 1 then "," else "") + jsonVector(values[i, :], significantDigits);
  end for;
  Modelica.Utilities.Streams.print(text + "]}", fileName);
end writeSliceField;
