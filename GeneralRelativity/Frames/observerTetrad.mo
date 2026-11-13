within GeneralRelativity.Frames;
function observerTetrad "Orthonormal tetrad whose timelike leg is a given observer four-velocity; spatial legs by Gram-Schmidt of d_1, d_2, d_3"
  input Real g[4, 4] "Covariant metric";
  input Real u[4] "Contravariant observer four-velocity (timelike, future-pointing; normalized here)";
  output Real e[4, 4] "e[a, :] = contravariant components of frame vector a; e[1, :] = u / |u|";
protected
  Real norm;
  Real v[4];
algorithm
  norm := -(u * (g * u));
  assert(norm > 0.0, "observerTetrad: the observer four-velocity is not timelike");
  e[1, :] := u / sqrt(norm);
  for i in 2:4 loop
    v := zeros(4);
    v[i] := 1.0;
    // Modified Gram-Schmidt against the legs found so far; eta_11 = -1 flips the projection sign.
    for j in 1:i - 1 loop
      v := v - (if j == 1 then -1.0 else 1.0) * (v * (g * e[j, :])) * e[j, :];
    end for;
    norm := v * (g * v);
    assert(norm > 0.0, "observerTetrad: coordinate basis vectors are degenerate with the observer");
    e[i, :] := v / sqrt(norm);
  end for;
  annotation(Documentation(info="<html>
<p>Satisfies <code>e g e<sup>T</sup> = diag(&minus;1, 1, 1, 1)</code>. Leg <code>e<sub>1</sub></code> lies in the
plane of <code>u</code> and <code>&part;<sub>1</sub></code> (the radial direction in spherical charts),
<code>e<sub>2</sub></code> adds <code>&part;<sub>2</sub></code> and <code>e<sub>3</sub></code> completes the frame.
It fails only when some combination of <code>&part;<sub>1</sub>, &part;<sub>2</sub>, &part;<sub>3</sub></code>
is parallel to <code>u</code> (an observer with <code>u<sup>0</sup> = 0</code>).</p>
</html>"));
end observerTetrad;
