within GeneralRelativity.Examples;
model SchwarzschildOrbitScene "Schwarzschild and Newton-Cartan orbits with the Three.js scene stream: Flamm's paraboloid coloured by curvature, horizon, photon sphere, radial profile and both animated bodies"
  extends GeneralRelativity.Examples.SchwarzschildNewtonCartanOrbit;
  parameter String sceneFile = "schwarzschild_scene.ndjson" "Scene file name for Visualization/threejs/index.html";
  parameter String sceneDirectory = "" "Output folder; empty uses GR_SCENE_DIR, else ~/Documents/Visualization, else the working directory";
  parameter Real sampleInterval = 1.0 "Proper time between scene samples";
  parameter Real rSceneMax = 30.0 * M "Outer radius of the embedding surface and profile";
initial algorithm
  GeneralRelativity.Visualization.beginScene(GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory), "Schwarzschild vs Newton-Cartan orbit", "spherical", rs, "M", 1.0);
  GeneralRelativity.Visualization.writeSphere(GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory), "horizon", rs, zeros(3), "#000000");
  GeneralRelativity.Visualization.writeSphere(GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory), "photon sphere", 1.5 * rs, zeros(3), "#5fb3ff");
  GeneralRelativity.Visualization.writeEmbeddingSurface(GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory), metric, p, step, 1.05 * rs, rSceneMax, 40, 64);
  GeneralRelativity.Visualization.writeRadialProfile(GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory), metric, p, step, 1.05 * rs, rSceneMax, 64);
  GeneralRelativity.Visualization.writeTrack(GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory), "schwarzschild", "#ffb454",
    {"Kretschmann", "E", "L", "separation"}, "Schwarzschild geodesic (general relativity)");
  GeneralRelativity.Visualization.writeTrack(GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory), "newton", "#7ee0a1",
    {"E Newton", "L Newton"}, "Newton-Cartan geodesic (Newtonian gravity), same initial position and dx/dt");
algorithm
  when sample(0.0, sampleInterval) then
    GeneralRelativity.Visualization.writeSample(GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory), time, x,
      {k.kretschmann, energy, angularMomentum, separation}, "schwarzschild");
    GeneralRelativity.Visualization.writeSample(GeneralRelativity.Visualization.scenePath(sceneFile, sceneDirectory), newton.x[1], newton.x,
      {newtonEnergy, newtonAngularMomentum}, "newton");
  end when;
  annotation(Documentation(info="<html>
<p>Adds the newline-delimited JSON scene stream to
<a href=\"modelica://GeneralRelativity.Examples.SchwarzschildNewtonCartanOrbit\">SchwarzschildNewtonCartanOrbit</a>.
Two tracks are written: <code>schwarzschild</code> (the relativistic geodesic, sampled
in proper time) and <code>newton</code> (the Newton-Cartan geodesic). Both samples
carry the same coordinate time <code>x[1] = ct</code>, and the viewer animates the
bodies together against it, showing the relativistic orbit precess away from the
closed Newtonian ellipse.</p>
<p>The file goes to <code>~/Documents/Visualization</code> by default; see
<a href=\"modelica://GeneralRelativity.Visualization.scenePath\">scenePath</a>
for how to choose another folder. Open it with
<code>GeneralRelativity/Visualization/threejs/index.html</code>. The writers use
<code>Modelica.Utilities.Streams</code>, so this model requires an OpenModelica build whose
<code>ModelicaExternalC</code> library links; the physics model does not.</p>
</html>"));
end SchwarzschildOrbitScene;
