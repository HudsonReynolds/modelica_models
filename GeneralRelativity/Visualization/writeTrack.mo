within GeneralRelativity.Visualization;
impure function writeTrack "Declare a trajectory track: its colour, label and the names of the scalars its samples carry"
  input String fileName;
  input String name "Track name used by writeSample";
  input String color = "#ffb454" "CSS colour of the body and its path";
  input String scalarNames[:] = fill("", 0) "Labels for the scalars of each sample, in order";
  input String description = "" "Legend text";
protected
  String names;
algorithm
  names := "[";
  for i in 1:size(scalarNames, 1) loop
    names := names + (if i > 1 then "," else "") + "\"" + scalarNames[i] + "\"";
  end for;
  names := names + "]";
  Modelica.Utilities.Streams.print("{\"type\":\"track\",\"name\":\"" + name + "\",\"color\":\"" + color
    + "\",\"scalarNames\":" + names + ",\"description\":\"" + description + "\"}", fileName);
end writeTrack;
