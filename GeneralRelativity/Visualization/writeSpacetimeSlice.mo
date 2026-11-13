within GeneralRelativity.Visualization;
impure function writeSpacetimeSlice "Append the equatorial slice of a metric sampled on a square grid: lapse, river, frame dragging and light cones"
  input String fileName;
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric function";
  input Real p[:] "Metric parameters";
  input Boolean spherical = true "true for {ct, r, theta, phi}, false for {ct, x, y, z}";
  input Real extent = 30.0 "Half-width of the square grid, in chart length units";
  input Integer n = 81 "Grid nodes per side";
  input Real rMin = 0.0 "Nodes closer to the origin are marked invalid without evaluating the metric";
  input Real ct = 0.0 "Coordinate time x[1] of the slice";
  input String name = "equator" "Slice name";
  input Integer significantDigits = 9 "Digits written per value";
protected
  Real coordinate[n];
  Real valid[n, n];
  Real lapse[n, n];
  Real timeLag[n, n];
  Real river[n, n];
  Real drag[n, n];
  Real coneValid[n, n];
  Real coneOut[n, n];
  Real coneIn[n, n];
  Real conePro[n, n];
  Real coneRetro[n, n];
  Real radius;
  Real x[4];
  SpacetimePoint s;
algorithm
  for i in 1:n loop
    coordinate[i] := -extent + 2.0 * extent * (i - 1) / (n - 1);
  end for;
  for i in 1:n loop
    for j in 1:n loop
      radius := sqrt(coordinate[i]^2 + coordinate[j]^2);
      if radius <= max(rMin, 1.0e-9 * extent) then
        s := SpacetimePoint(valid = 0.0, lapse = 0.0, timeLag = 1.0, riverSpeed = 0.0, dragSpeed = 0.0,
          coneValid = 0.0, coneOutward = 0.0, coneInward = 0.0, conePrograde = 0.0, coneRetrograde = 0.0);
      else
        x := if spherical then {ct, radius, GeneralRelativity.Constants.pi / 2.0, atan2(coordinate[j], coordinate[i])}
          else {ct, coordinate[i], coordinate[j], 0.0};
        s := spacetimePointAt(metric, x, p, spherical);
      end if;
      valid[i, j] := s.valid;
      lapse[i, j] := s.lapse;
      timeLag[i, j] := s.timeLag;
      river[i, j] := s.riverSpeed;
      drag[i, j] := s.dragSpeed;
      coneValid[i, j] := s.coneValid;
      coneOut[i, j] := s.coneOutward;
      coneIn[i, j] := s.coneInward;
      conePro[i, j] := s.conePrograde;
      coneRetro[i, j] := s.coneRetrograde;
    end for;
  end for;
  Modelica.Utilities.Streams.print("{\"type\":\"slice\",\"name\":\"" + name + "\",\"plane\":\"xy\",\"chart\":\""
    + (if spherical then "spherical" else "cartesian") + "\",\"ct\":" + jsonNumber(ct) + ",\"rMin\":" + jsonNumber(rMin)
    + ",\"u\":" + jsonVector(coordinate) + ",\"v\":" + jsonVector(coordinate) + "}", fileName);
  writeSliceField(fileName, name, "valid", valid, 1);
  writeSliceField(fileName, name, "lapse", lapse, significantDigits);
  writeSliceField(fileName, name, "timeLag", timeLag, significantDigits);
  writeSliceField(fileName, name, "river", river, significantDigits);
  writeSliceField(fileName, name, "drag", drag, significantDigits);
  writeSliceField(fileName, name, "coneValid", coneValid, 1);
  writeSliceField(fileName, name, "coneOut", coneOut, significantDigits);
  writeSliceField(fileName, name, "coneIn", coneIn, significantDigits);
  writeSliceField(fileName, name, "conePro", conePro, significantDigits);
  writeSliceField(fileName, name, "coneRetro", coneRetro, significantDigits);
  annotation(Documentation(info="<html>
<p>Writes a <code>slice</code> record with the grid coordinates, then one <code>sliceField</code>
record per quantity of <a href=\"modelica://GeneralRelativity.Visualization.SpacetimePoint\">SpacetimePoint</a>,
<code>values[i][j]</code> at <code>(u[i], v[j])</code>. The plane is the equator <code>&theta; = &pi;/2</code> of a
spherical chart or <code>z = 0</code> of a Cartesian one: one spatial dimension is dropped, which is what
the river, clock and spacetime-diagram styles of the viewer draw. Add further fields with
<a href=\"modelica://GeneralRelativity.Visualization.writeSliceField\">writeSliceField</a>.</p>
<p>Time-dependent metrics (the cislunar metric with its moving Moon) are sampled at <code>ct</code> only.</p>
</html>"));
end writeSpacetimeSlice;
