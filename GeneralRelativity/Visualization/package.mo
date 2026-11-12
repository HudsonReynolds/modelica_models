within GeneralRelativity;
package Visualization "Embedding surfaces, curvature fields and a newline-delimited JSON scene stream for the Three.js viewer"
  annotation(Documentation(info="<html>
<h4>Scene stream</h4>
<p>Models write one JSON object per line to a text file through
<code>Modelica.Utilities.Streams.print</code>. The viewer in
<code>GeneralRelativity/Visualization/threejs/index.html</code> opens that file
in a browser without a server and renders:</p>
<ul>
<li><code>scene</code>: name, chart (<code>spherical</code> or <code>cartesian</code>),
length unit and a reference radius (<a href=\"modelica://GeneralRelativity.Visualization.beginScene\">beginScene</a>);</li>
<li><code>surfaceRow</code>: one radius of an embedding surface with its Kretschmann
colour value (<a href=\"modelica://GeneralRelativity.Visualization.writeEmbeddingSurface\">writeEmbeddingSurface</a>);</li>
<li><code>profile</code>: <code>g<sub>tt</sub></code>, <code>g<sub>rr</sub></code>, Ricci scalar and
Kretschmann against radius (<a href=\"modelica://GeneralRelativity.Visualization.writeRadialProfile\">writeRadialProfile</a>);</li>
<li><code>field</code>: a scalar sampled on a plane, for tidal or curvature maps
(<a href=\"modelica://GeneralRelativity.Visualization.writeFieldSlice\">writeFieldSlice</a>);</li>
<li><code>sphere</code>: a reference sphere such as a horizon or planet
(<a href=\"modelica://GeneralRelativity.Visualization.writeSphere\">writeSphere</a>);</li>
<li><code>sample</code>: one point of a trajectory with scalars attached
(<a href=\"modelica://GeneralRelativity.Visualization.writeSample\">writeSample</a>).</li>
</ul>
<h4>Embedding surfaces</h4>
<p>For a static spherically symmetric metric with <code>g<sub>&theta;&theta;</sub> = r&sup2;</code>
the equatorial slice <code>t = const, &theta; = &pi;/2</code> embeds in Euclidean
3-space as a surface of revolution with <code>dz/dr = sqrt(g<sub>rr</sub> - 1)</code>.
For Schwarzschild this is Flamm's paraboloid <code>z = 2 sqrt(r<sub>s</sub>(r - r<sub>s</sub>))</code>
(<a href=\"modelica://GeneralRelativity.Visualization.flammParaboloid\">flammParaboloid</a>).
The surface is coloured by the Kretschmann scalar so curvature is visible
even where the embedding is nearly flat.</p>
<p>The writers are impure and may only be called from initial sections,
when-clauses or other impure functions.</p>
</html>"));
end Visualization;
