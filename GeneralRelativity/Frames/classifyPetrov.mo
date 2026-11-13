within GeneralRelativity.Frames;
function classifyPetrov "Petrov type of a Weyl tensor at an event"
  input Real g[4, 4] "Covariant metric (g^00 < 0 for the Eulerian frame)";
  input Real C[4, 4, 4, 4] "Weyl tensor C_abcd";
  input Real curvatureScale "Reference curvature magnitude; a Weyl tensor below tolerance * curvatureScale is type O";
  input Real tolerance = 1.0e-3 "Relative tolerance of the degeneracy tests";
  output PetrovType petrov "Petrov type";
protected
  Real E[3, 3];
  Real B[3, 3];
algorithm
  (E, B) := weylElectricMagnetic(g, C);
  petrov := classifyPetrovMatrix(E, B, curvatureScale, tolerance);
  annotation(Documentation(info="<html>
<p>Builds <code>Q = E + iB</code> with
<a href=\"modelica://GeneralRelativity.Frames.weylElectricMagnetic\">weylElectricMagnetic</a> and classifies it with
<a href=\"modelica://GeneralRelativity.Frames.classifyPetrovMatrix\">classifyPetrovMatrix</a>. The Petrov type does not
depend on the observer, so the Eulerian frame is as good as any. A natural <code>curvatureScale</code> is
the frame norm of the Riemann tensor, as used by
<a href=\"modelica://GeneralRelativity.Frames.classifyPetrovCurvature\">classifyPetrovCurvature</a>; pass
<code>0</code> to report <b>O</b> only for an exactly vanishing Weyl tensor.</p>
</html>"));
end classifyPetrov;
