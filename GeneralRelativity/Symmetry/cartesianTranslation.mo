within GeneralRelativity.Symmetry;
function cartesianTranslation "Constant vector field p = {p0, px, py, pz} in a Cartesian chart"
  extends GeneralRelativity.Interfaces.VectorFieldFunction;
algorithm
  xi := p[1:4];
  annotation(Documentation(info="<html>
<p>Spacetime translations are Killing vectors of Minkowski space in
<a href=\"modelica://GeneralRelativity.Metrics.minkowski\">Cartesian coordinates</a>.</p>
</html>"));
end cartesianTranslation;
