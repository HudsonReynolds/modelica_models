within GeneralRelativity.Visualization;
impure function writeSphere "Append a reference sphere (horizon, photon sphere, planet) to the scene"
  input String fileName;
  input String name;
  input Real radius;
  input Real centre[3] = {0.0, 0.0, 0.0};
  input String color = "#888888" "CSS colour";
algorithm
  Modelica.Utilities.Streams.print("{\"type\":\"sphere\",\"name\":\"" + name + "\",\"radius\":" + jsonNumber(radius)
    + ",\"centre\":" + jsonVector(centre) + ",\"color\":\"" + color + "\"}", fileName);
end writeSphere;
