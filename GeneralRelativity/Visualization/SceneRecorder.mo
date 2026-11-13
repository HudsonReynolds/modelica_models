within GeneralRelativity.Visualization;
model SceneRecorder "Reusable scene writer: one component gives a simulation every viewer style (embedding, river, clocks, spacetime diagram, ray tracer)"
  replaceable function metric = GeneralRelativity.Metrics.minkowski
    constrainedby GeneralRelativity.Interfaces.MetricFunction "Metric of the simulation" annotation(choicesAllMatching = true);
  parameter Real p[:] = {0.0} "Metric parameters";
  parameter Real step[4] = {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3} "Finite-difference step for the embedding and profile curvature";
  parameter Boolean spherical = true "Chart of the metric and samples: true for {ct, r, theta, phi}, false for {ct, x, y, z}";
  parameter String metricName = "" "Analytic family for the ray tracer (schwarzschild, kerr, reissnerNordstrom, kerrNewman, minkowski); empty if none applies";

  parameter String sceneFile = "scene.ndjson" "Scene file name";
  parameter String sceneDirectory = "" "Output folder; empty uses GR_SCENE_DIR, else ~/Documents/Visualization, else the working directory";
  parameter String title = "scene" "Scene title shown by the viewer";
  parameter Real referenceRadius = 1.0 "Length that frames the view, for example the Schwarzschild radius";
  parameter String lengthUnit = "M" "Label for lengths";
  parameter Real c = 1.0 "Speed of light in chart units";
  parameter Real rInner = 1.0 "Innermost radius drawn, outside any horizon or body surface";
  parameter Real rOuter = 30.0 "Outer radius of the embedding, profile and slice";

  parameter Boolean writeEmbedding = spherical "Embedding surface of the equatorial plane (spherical charts)" annotation(Dialog(group = "Content"));
  parameter Boolean writeProfile = spherical "Radial profile of g_tt, g_rr, R and K (spherical charts)" annotation(Dialog(group = "Content"));
  parameter Boolean writeSlice = true "Equatorial slice for the river, clock and spacetime styles" annotation(Dialog(group = "Content"));
  parameter Integer sliceResolution = 81 "Slice grid nodes per side" annotation(Dialog(group = "Content"));
  parameter Real sliceInner = rInner "Innermost slice radius; below a horizon the river and its light flashes are still drawn, the lapse is marked invalid" annotation(Dialog(group = "Content"));
  parameter String sphereNames[:] = fill("", 0) "Reference spheres, e.g. {\"horizon\", \"photon sphere\"}" annotation(Dialog(group = "Content"));
  parameter Real sphereRadii[size(sphereNames, 1)] = fill(1.0, size(sphereNames, 1)) "Sphere radii" annotation(Dialog(group = "Content"));
  parameter String sphereColors[size(sphereNames, 1)] = fill("#888888", size(sphereNames, 1)) "CSS colours" annotation(Dialog(group = "Content"));

  parameter Integer nTracks = 1 "Number of bodies" annotation(Dialog(group = "Tracks"));
  parameter String trackNames[nTracks] = {"body" + String(i) for i in 1:nTracks} "Body names" annotation(Dialog(group = "Tracks"));
  parameter String trackColors[nTracks] = fill("#ffb454", nTracks) "Body colours" annotation(Dialog(group = "Tracks"));
  parameter String trackDescriptions[nTracks] = fill("", nTracks) "Legend text" annotation(Dialog(group = "Tracks"));
  parameter Integer nScalars = 0 "Scalars attached to every sample" annotation(Dialog(group = "Tracks"));
  parameter String scalarNames[nTracks, nScalars] = fill("", nTracks, nScalars) "Scalar labels per body; empty labels are hidden" annotation(Dialog(group = "Tracks"));
  parameter Real sampleInterval = 1.0 "Modelica time between samples" annotation(Dialog(group = "Tracks"));

  parameter String style = "embedding" "Style the viewer opens with: embedding, river, clocks, spacetime or raytrace" annotation(Dialog(group = "View"),
    choices(choice = "\"embedding\"", choice = "\"river\"", choice = "\"clocks\"", choice = "\"spacetime\"", choice = "\"raytrace\""));
  parameter String optionNames[:] = fill("", 0) "Style options to preset" annotation(Dialog(group = "View"));
  parameter Real optionValues[size(optionNames, 1)] = fill(0.0, size(optionNames, 1)) "Their values" annotation(Dialog(group = "View"));

  input Real trackX[nTracks, 4] "Coordinates of each body; x[1] = c t is the clock the viewer animates against";
  input Real trackLambda[nTracks] "Own clock of each body in time units: proper time for a massive body, absolute time t for a Newtonian one; the clocks style compares it with t = x[1]/c";
  input Real trackScalars[nTracks, nScalars] = zeros(nTracks, nScalars) "Scalars of each body";
initial algorithm
  beginScene(scenePath(sceneFile, sceneDirectory), title, if spherical then "spherical" else "cartesian", referenceRadius, lengthUnit, c);
  writeView(scenePath(sceneFile, sceneDirectory), style, optionNames, optionValues);
  if metricName <> "" then
    writeMetric(scenePath(sceneFile, sceneDirectory), metricName, p, if spherical then "spherical" else "cartesian");
  end if;
  for i in 1:size(sphereNames, 1) loop
    writeSphere(scenePath(sceneFile, sceneDirectory), sphereNames[i], sphereRadii[i], zeros(3), sphereColors[i]);
  end for;
  if writeEmbedding then
    writeEmbeddingSurface(scenePath(sceneFile, sceneDirectory), metric, p, step, rInner, rOuter, 40, 64);
  end if;
  if writeProfile then
    writeRadialProfile(scenePath(sceneFile, sceneDirectory), metric, p, step, rInner, rOuter, 64);
  end if;
  if writeSlice then
    writeSpacetimeSlice(scenePath(sceneFile, sceneDirectory), metric, p, spherical, rOuter, sliceResolution, sliceInner, trackX[1, 1]);
  end if;
  for i in 1:nTracks loop
    writeTrack(scenePath(sceneFile, sceneDirectory), trackNames[i], trackColors[i], scalarNames[i, :], trackDescriptions[i]);
  end for;
algorithm
  when sample(0.0, sampleInterval) then
    for i in 1:nTracks loop
      writeSample(scenePath(sceneFile, sceneDirectory), trackLambda[i], trackX[i, :], trackScalars[i, :], trackNames[i]);
    end for;
  end when;
  annotation(Documentation(info="<html>
<p>Drop one of these into any simulation to get a scene file for
<code>Visualization/threejs/index.html</code>; the same file can then be viewed in every style the viewer
offers, switched at run time:</p>
<table>
<tr><th>Style</th><th>Shows</th><th>Needs</th></tr>
<tr><td><code>embedding</code></td><td>3D view; the equatorial plane as an embedded surface coloured by curvature</td><td>nothing (surface: <code>writeEmbedding</code>)</td></tr>
<tr><td><code>river</code></td><td>River model: space flowing inward at <code>&radic;(1+g<sub>00</sub>)</code>, frame dragging, light flashes carried by the flow</td><td><code>writeSlice</code></td></tr>
<tr><td><code>clocks</code></td><td>Time-contracting grid: a sheet that sinks by the accumulated time lag <code>(1&minus;&alpha;)t</code>, with ticking clocks and body clocks</td><td><code>writeSlice</code></td></tr>
<tr><td><code>spacetime</code></td><td>2+1 spacetime diagram: worldlines against coordinate time, light cones that close toward the horizon</td><td>tracks; cones need <code>writeSlice</code></td></tr>
<tr><td><code>raytrace</code></td><td>GPU ray tracer of null geodesics with the bodies, an accretion disk and a lensed sky</td><td><code>metricName</code></td></tr>
</table>
<p>Example, inside a model that integrates two bodies:</p>
<pre>GeneralRelativity.Visualization.SceneRecorder scene(
  redeclare function metric = GeneralRelativity.Metrics.schwarzschild,
  p = {rs}, metricName = \"schwarzschild\", sceneFile = \"orbit.ndjson\",
  referenceRadius = rs, rInner = 1.05*rs, rOuter = 30*M,
  sphereNames = {\"horizon\"}, sphereRadii = {rs}, sphereColors = {\"#000000\"},
  nTracks = 2, trackNames = {\"gr\", \"newton\"},
  trackX = {x, newton.x}, trackLambda = {time, newton.x[1]},
  style = \"river\");</pre>
<p>Bind <code>trackScalars</code> to plain variables. A field of a record returned by a function, such as
<code>k.kretschmann</code> of <code>Metrics.curvatureAt</code>, should first be held in a <code>discrete</code>
variable updated in <code>when sample(0, sampleInterval)</code>: bound continuously, OpenModelica 1.28
re-evaluates the function many times per step or fails to scalarize it.</p>
<p>All tracks should share <code>x[1] = ct</code>; the viewer animates them against it. The slice is
sampled at the first track's <code>x[1]</code> at the start, so time-dependent metrics are shown as a
snapshot.</p>
</html>"));
end SceneRecorder;
