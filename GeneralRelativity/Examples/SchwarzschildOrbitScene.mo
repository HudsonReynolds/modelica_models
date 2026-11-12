within GeneralRelativity.Examples;
model SchwarzschildOrbitScene "Schwarzschild orbit with the Three.js scene stream: Flamm's paraboloid coloured by curvature, horizon, photon sphere, radial profile and the sampled orbit"
  extends GeneralRelativity.Examples.SchwarzschildOrbit;
  parameter String sceneFile = "schwarzschild_scene.ndjson" "Scene stream for Visualization/threejs/index.html, written in the working directory";
  parameter Real sampleInterval = 2.0 "Proper time between scene samples";
  parameter Real rSceneMax = 30.0 * M "Outer radius of the embedding surface and profile";
initial algorithm
  GeneralRelativity.Visualization.beginScene(sceneFile, "Schwarzschild orbit", "spherical", rs, "M", 1.0);
  GeneralRelativity.Visualization.writeSphere(sceneFile, "horizon", rs, zeros(3), "#000000");
  GeneralRelativity.Visualization.writeSphere(sceneFile, "photon sphere", 1.5 * rs, zeros(3), "#5fb3ff");
  GeneralRelativity.Visualization.writeEmbeddingSurface(sceneFile, metric, p, step, 1.05 * rs, rSceneMax, 40, 64);
  GeneralRelativity.Visualization.writeRadialProfile(sceneFile, metric, p, step, 1.05 * rs, rSceneMax, 64);
algorithm
  when sample(0.0, sampleInterval) then
    GeneralRelativity.Visualization.writeSample(sceneFile, time, x, {k.kretschmann, energy, angularMomentum});
  end when;
  annotation(Documentation(info="<html>
<p>Adds the newline-delimited JSON scene stream to
<a href=\"modelica://GeneralRelativity.Examples.SchwarzschildOrbit\">SchwarzschildOrbit</a>. Open the file
with <code>GeneralRelativity/Visualization/threejs/index.html</code>. The writers use
<code>Modelica.Utilities.Streams</code>, so this model requires an OpenModelica build whose
<code>ModelicaExternalC</code> library links; the physics model does not.</p>
</html>"));
end SchwarzschildOrbitScene;
