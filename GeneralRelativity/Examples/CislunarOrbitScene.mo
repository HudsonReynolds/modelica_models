within GeneralRelativity.Examples;
model CislunarOrbitScene "Cislunar orbit with the Three.js scene stream: Earth, the Moon at t = 0 and the sampled trajectory"
  extends GeneralRelativity.Examples.CislunarOrbit;
  parameter String sceneFile = "cislunar_scene.ndjson" "Scene stream for Visualization/threejs/index.html, written in the working directory";
  parameter Real sampleInterval = 1800.0 "Seconds between scene samples";
initial algorithm
  GeneralRelativity.Visualization.beginScene(sceneFile, "Cislunar orbit", "cartesian", GeneralRelativity.Constants.R_earth, "m", c);
  GeneralRelativity.Visualization.writeSphere(sceneFile, "Earth", GeneralRelativity.Constants.R_earth, zeros(3), "#3b7dd8");
  GeneralRelativity.Visualization.writeSphere(sceneFile, "Moon at t = 0", 1737.4e3,
    GeneralRelativity.Perturbation.Potentials.circularBodyPosition(0.0, moonDistance, moonRate, moonPhase0), "#bbbbbb");
algorithm
  when sample(0.0, sampleInterval) then
    GeneralRelativity.Visualization.writeSample(sceneFile, time, {c * time, x[1], x[2], x[3]}, {jacobi - jacobi0, moonRange});
  end when;
  annotation(Documentation(info="<html>
<p>Adds the newline-delimited JSON scene stream to
<a href=\"modelica://GeneralRelativity.Examples.CislunarOrbit\">CislunarOrbit</a>. Open the file
with <code>GeneralRelativity/Visualization/threejs/index.html</code>. The writers use
<code>Modelica.Utilities.Streams</code>, so this model requires an OpenModelica build whose
<code>ModelicaExternalC</code> library links; the physics model does not.</p>
</html>"));
end CislunarOrbitScene;
