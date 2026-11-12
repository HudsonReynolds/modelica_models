within GeneralRelativity.StressEnergy;
function cosmologicalConstant "Vacuum energy T_ab = -Lambda c^4 / (8 pi G) g_ab"
  input Real g[4, 4] "Covariant metric";
  input Real Lambda(unit = "1/m2") "Cosmological constant";
  input Real gravitationalConstant(unit = "m3/(kg.s2)") = 1.0 "Newton's constant";
  input Real c(unit = "m/s") = 1.0 "Speed of light";
  output Real T[4, 4] "T_ab";
algorithm
  T := -Lambda * c^4 / (8.0 * GeneralRelativity.Constants.pi * gravitationalConstant) * g;
  annotation(Documentation(info="<html>
<p>Moving <code>&Lambda;g<sub>ab</sub></code> to the right-hand side of the field
equations gives a perfect fluid with <code>p = -&rho;c&sup2;</code>.</p>
</html>"));
end cosmologicalConstant;
