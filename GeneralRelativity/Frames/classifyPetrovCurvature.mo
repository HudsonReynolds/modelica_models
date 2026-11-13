within GeneralRelativity.Frames;
function classifyPetrovCurvature "Petrov type of the Weyl tensor of a Curvature record, scaled by the Riemann tensor"
  input GeneralRelativity.Tensors.Curvature k "Curvature at one event, e.g. from Metrics.curvatureAt";
  input Real tolerance = 1.0e-3 "Relative tolerance of the degeneracy tests";
  output PetrovType petrov "Petrov type";
protected
  Real F[4, 4, 4, 4] "Riemann tensor in the Eulerian orthonormal frame";
  Real scale;
algorithm
  F := frameComponents4(k.riemann, orthonormalTetrad(k.g));
  // sqrt(sum F^2 / 8) equals |E + iB| when the Riemann tensor is pure Weyl.
  scale := sqrt(sum(F[a, b, c, d]^2 for a in 1:4, b in 1:4, c in 1:4, d in 1:4) / 8.0);
  petrov := classifyPetrov(k.g, GeneralRelativity.Tensors.weyl(k.g, k.riemann, k.ricci, k.ricciScalar), scale, tolerance);
  annotation(Documentation(info="<html>
<p>The type-O threshold is <code>tolerance</code> times the Euclidean norm of the orthonormal-frame Riemann
components divided by <code>&radic;8</code>, which equals <code>|E + iB|</code> in vacuum. A conformally flat
metric with nonzero Ricci curvature (such as FLRW) is therefore reported as <b>O</b> despite finite-difference
noise in its Weyl tensor.</p>
</html>"));
end classifyPetrovCurvature;
