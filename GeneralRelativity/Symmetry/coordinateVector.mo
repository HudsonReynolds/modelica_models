within GeneralRelativity.Symmetry;
function coordinateVector "Coordinate basis vector d_k; p = {k} with k = 1..4"
  extends GeneralRelativity.Interfaces.VectorFieldFunction;
algorithm
  xi := zeros(4);
  xi[integer(p[1] + 0.5)] := 1.0;
  annotation(Documentation(info="<html>
<p><code>p = {1}</code> is <code>&part;<sub>0</sub> = &part;/&part;(ct)</code>,
<code>p = {2}</code> is <code>&part;<sub>r</sub></code> in a spherical chart, and so on.
A coordinate vector is Killing exactly when the metric does not depend on that
coordinate.</p>
</html>"));
end coordinateVector;
