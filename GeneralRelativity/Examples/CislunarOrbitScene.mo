within GeneralRelativity.Examples;
model CislunarOrbitScene "Cislunar orbit with the Three.js scene stream: Earth, the Moon at t = 0 and the sampled trajectory"
  extends GeneralRelativity.Examples.CislunarOrbit;
  parameter String sceneFile = "cislunar_scene.ndjson" "Scene file name for Visualization/threejs/index.html";
  parameter String sceneDirectory = "" "Output folder; empty uses GR_SCENE_DIR, else ~/Documents/Visualization, else the working directory";
  parameter Real sampleInterval = 1800.0 "Seconds between scene samples";
initial algorithm
  GeneralRelativity.Visualization.beginScene(GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory), "Cislunar orbit", "cartesian", GeneralRelativity.Constants.R_earth, "m", c);
  GeneralRelativity.Visualization.writeSphere(GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory), "Earth", GeneralRelativity.Constants.R_earth, zeros(3), "#3b7dd8");
  GeneralRelativity.Visualization.writeSphere(GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory), "Moon at t = 0", 1737.4e3,
    GeneralRelativity.Perturbation.Potentials.circularBodyPosition(0.0, moonDistance, moonRate, moonPhase0), "#bbbbbb");
  GeneralRelativity.Visualization.writeView(GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory), "embedding");
  // Earth-Moon metric at t = 0 (the Moon moves, so this is a snapshot); the viewer exaggerates the weak field.
  GeneralRelativity.Visualization.writeSpacetimeSlice(GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory),
    GeneralRelativity.Metrics.cislunar, metricParameters, false, 1.15 * moonDistance, 81, GeneralRelativity.Constants.R_earth, 0.0);
algorithm
  when sample(0.0, sampleInterval) then
    GeneralRelativity.Visualization.writeSample(GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory), time, {c * time, x[1], x[2], x[3]}, {jacobi - jacobi0, moonRange});
  end when;
  annotation(Documentation(info="<html>
<p>Adds the newline-delimited JSON scene stream to
<a href=\"modelica://GeneralRelativity.Examples.CislunarOrbit\">CislunarOrbit</a>. The file goes to <code>~/Documents/Visualization</code> by default; see
<a href=\"modelica://GeneralRelativity.Visualization.scenePath\">scenePath</a>
for how to choose another folder. Open it with
<code>GeneralRelativity/Visualization/threejs/index.html</code>. The writers use
<code>Modelica.Utilities.Streams</code>, so this model requires an OpenModelica build whose
<code>ModelicaExternalC</code> library links; the physics model does not.</p>
</html>"));
end CislunarOrbitScene;
