within GeneralRelativity.Examples;
model SchwarzschildOrbitScene "Schwarzschild and Newton-Cartan orbits recorded for the Three.js viewer, viewable as embedding, river, clocks, spacetime diagram or ray tracer"
  extends GeneralRelativity.Examples.SchwarzschildNewtonCartanOrbit;
  parameter String sceneFile = "schwarzschild_scene.ndjson" "Scene file name for Visualization/threejs/index.html";
  parameter String sceneDirectory = "" "Output folder; empty uses GR_SCENE_DIR, else ~/Documents/Visualization, else the working directory";
  parameter String style = "river" "Style the viewer opens with; every style can be chosen in the viewer"
    annotation(choices(choice = "\"embedding\"", choice = "\"river\"", choice = "\"clocks\"", choice = "\"spacetime\"", choice = "\"raytrace\""));
  parameter Real sampleInterval = 1.0 "Proper time between scene samples";
  parameter Real rSceneMax = 30.0 * M "Outer radius of the embedding surface, profile and slice";
  discrete Real kretschmannSampled(start = 0.0, fixed = true) "k.kretschmann held at the scene samples";
  GeneralRelativity.Visualization.SceneRecorder scene(
    redeclare function metric = GeneralRelativity.Metrics.schwarzschild,
    p = {rs},
    step = step,
    spherical = true,
    metricName = "schwarzschild",
    sceneFile = sceneFile,
    sceneDirectory = sceneDirectory,
    title = "Schwarzschild vs Newton-Cartan orbit",
    referenceRadius = rs,
    lengthUnit = "M",
    c = c,
    rInner = 1.05 * rs,
    rOuter = rSceneMax,
    sliceInner = 0.25 * rs,
    sphereNames = {"horizon", "photon sphere"},
    sphereRadii = {rs, 1.5 * rs},
    sphereColors = {"#000000", "#5fb3ff"},
    nTracks = 2,
    trackNames = {"schwarzschild", "newton"},
    trackColors = {"#ffb454", "#7ee0a1"},
    trackDescriptions = {"Schwarzschild geodesic (general relativity)", "Newton-Cartan geodesic (Newtonian gravity), same initial position and dx/dt"},
    nScalars = 4,
    scalarNames = {{"Kretschmann", "E", "L", "separation"}, {"E Newton", "L Newton", "", ""}},
    sampleInterval = sampleInterval,
    style = style,
    trackX = {x, newton.x},
    trackLambda = {time, newton.t},
    trackScalars = {{kretschmannSampled, energy, angularMomentum, separation}, {newtonEnergy, newtonAngularMomentum, 0.0, 0.0}})
    "Writes the scene; both bodies share x[1] = ct";
equation
  // A field of the record returned by curvatureAt, bound continuously into another component, makes
  // OpenModelica re-evaluate the curvature many times per step (or fail to scalarize it); holding it at the
  // sample instants, the only time the recorder reads it, avoids that.
  when sample(0.0, sampleInterval) then
    kretschmannSampled = k.kretschmann;
  end when;
  annotation(Documentation(info="<html>
<p>Adds a <a href=\"modelica://GeneralRelativity.Visualization.SceneRecorder\">SceneRecorder</a> to
<a href=\"modelica://GeneralRelativity.Examples.SchwarzschildNewtonCartanOrbit\">SchwarzschildNewtonCartanOrbit</a>.
Two tracks are written: <code>schwarzschild</code> (the relativistic geodesic, sampled in proper time)
and <code>newton</code> (the Newton-Cartan geodesic). Both carry the same coordinate time
<code>x[1] = ct</code>, and the viewer animates them together against it.</p>
<p>The scene contains everything every viewer style needs, so the same file can be shown as</p>
<ul>
<li><b>river</b> (the default): space flowing into the hole at <code>&radic;(r<sub>s</sub>/r)</code>, faster
than light inside the horizon, with light flashes swept along by the flow;</li>
<li><b>clocks</b>: a grid that sinks by the time each static clock has lost, with the two bodies'
own clocks;</li>
<li><b>spacetime</b>: the worldlines in a 2+1 diagram with light cones that close toward the horizon;</li>
<li><b>embedding</b>: Flamm's paraboloid coloured by curvature;</li>
<li><b>raytrace</b>: the hole and both bodies as a camera would see them, with lensing.</li>
</ul>
<p>The file goes to <code>~/Documents/Visualization</code> by default; see
<a href=\"modelica://GeneralRelativity.Visualization.scenePath\">scenePath</a>
for how to choose another folder. Open it with
<code>GeneralRelativity/Visualization/threejs/index.html</code>. The writers use
<code>Modelica.Utilities.Streams</code>, so this model requires an OpenModelica build whose
<code>ModelicaExternalC</code> library links; the physics model does not.</p>
</html>"));
end SchwarzschildOrbitScene;
