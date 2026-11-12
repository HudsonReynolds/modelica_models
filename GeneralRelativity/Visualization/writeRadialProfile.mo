within GeneralRelativity.Visualization;
impure function writeRadialProfile "Append g_tt, g_rr, Ricci scalar and Kretschmann scalar along the equatorial radius"
  input String fileName;
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric in {ct, r, theta, phi}";
  input Real p[:] "Metric parameters";
  input Real step[4] "Finite-difference step for curvature";
  input Real rMin;
  input Real rMax;
  input Integer n = 64;
protected
  Real r[n];
  Real gtt[n];
  Real grr[n];
  Real ricci[n];
  Real kretschmann[n];
  GeneralRelativity.Tensors.Curvature k;
algorithm
  for i in 1:n loop
    r[i] := rMin + (rMax - rMin) * (i - 1) / (n - 1);
    k := GeneralRelativity.Metrics.curvatureAt(metric, {0.0, r[i], GeneralRelativity.Constants.pi / 2.0, 0.0}, p, step);
    gtt[i] := k.g[1, 1];
    grr[i] := k.g[2, 2];
    ricci[i] := k.ricciScalar;
    kretschmann[i] := k.kretschmann;
  end for;
  Modelica.Utilities.Streams.print("{\"type\":\"profile\",\"r\":" + jsonVector(r) + ",\"gtt\":" + jsonVector(gtt)
    + ",\"grr\":" + jsonVector(grr) + ",\"ricciScalar\":" + jsonVector(ricci)
    + ",\"kretschmann\":" + jsonVector(kretschmann) + "}", fileName);
end writeRadialProfile;
