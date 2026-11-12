within GeneralRelativity.Examples;
model EarthJ2OrbitScene "Earth J2 orbit with the Three.js scene stream: Earth sphere, equatorial tidal-magnitude slice and the sampled orbit"
  extends GeneralRelativity.Examples.EarthJ2Orbit;
  parameter String sceneFile = "earth_scene.ndjson" "Scene stream for Visualization/threejs/index.html, written in the working directory";
  parameter Real sampleInterval = 60.0 "Seconds between scene samples";
  final parameter Real sliceCoordinate[40] = {-3.0 * R + 6.0 * R * (i - 1) / 39.0 for i in 1:40} "Even count so no sample sits on the singular origin";
initial algorithm
  GeneralRelativity.Visualization.beginScene(sceneFile, "Earth J2 orbit", "cartesian", R, "m", c);
  GeneralRelativity.Visualization.writeSphere(sceneFile, "Earth", R, zeros(3), "#3b7dd8");
  GeneralRelativity.Visualization.writeFieldSlice(sceneFile, "tidal magnitude |Phi_,ij|", "xy", 0.0,
    sliceCoordinate, sliceCoordinate,
    GeneralRelativity.Visualization.tidalMagnitudeSlice(GeneralRelativity.Perturbation.Potentials.oblateJ2Gradient, p, 1.0, sliceCoordinate, sliceCoordinate, 0.0, R));
algorithm
  when sample(0.0, sampleInterval) then
    GeneralRelativity.Visualization.writeSample(sceneFile, time, {c * time, x[1], x[2], x[3]},
      {relativeCorrection, raan, tidalLaplacian});
  end when;
  annotation(Documentation(info="<html>
<p>Adds the newline-delimited JSON scene stream to
<a href=\"modelica://GeneralRelativity.Examples.EarthJ2Orbit\">EarthJ2Orbit</a>. Open the file
with <code>GeneralRelativity/Visualization/threejs/index.html</code>. The writers use
<code>Modelica.Utilities.Streams</code>, so this model requires an OpenModelica build whose
<code>ModelicaExternalC</code> library links; the physics model does not.</p>
</html>"));
end EarthJ2OrbitScene;
