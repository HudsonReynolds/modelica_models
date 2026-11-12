within GeneralRelativity.Visualization;
impure function beginScene "Start a scene file: delete any previous file and write the header line"
  input String fileName "Scene file, newline-delimited JSON";
  input String name "Scene title";
  input String chart "\"spherical\" for {ct, r, theta, phi} samples or \"cartesian\" for {ct, x, y, z}";
  input Real referenceRadius "Radius that frames the view, for example the Schwarzschild radius or a planet radius";
  input String lengthUnit = "M" "Label for lengths";
  input Real c = 1.0 "Speed of light in chart units, used to convert x[1] to time";
algorithm
  Modelica.Utilities.Files.removeFile(fileName);
  Modelica.Utilities.Streams.print("{\"type\":\"scene\",\"name\":\"" + name
    + "\",\"chart\":\"" + chart + "\",\"referenceRadius\":" + jsonNumber(referenceRadius)
    + ",\"lengthUnit\":\"" + lengthUnit + "\",\"c\":" + jsonNumber(c) + "}", fileName);
end beginScene;
