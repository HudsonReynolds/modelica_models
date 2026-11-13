within GeneralRelativity.Frames;
function orthonormalTetrad "Orthonormal tetrad of the Eulerian (ZAMO) observer normal to t = const, spatial legs by Gram-Schmidt of d_1, d_2, d_3"
  input Real g[4, 4] "Covariant metric";
  output Real e[4, 4] "e[a, :] = contravariant components of frame vector a; e[1, :] = unit normal";
protected
  Real gInv[4, 4];
  Boolean accepted;
algorithm
  (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(g);
  assert(accepted, "orthonormalTetrad: the metric is singular");
  assert(gInv[1, 1] < 0.0, "orthonormalTetrad: the slice t = const is not spacelike (g^00 >= 0)");
  e := observerTetrad(g, -gInv[:, 1] / sqrt(-gInv[1, 1]));
  annotation(Documentation(info="<html>
<p>The timelike leg is the future unit normal <code>n<sup>&mu;</sup> = &minus;g<sup>&mu;0</sup>/&radic;(&minus;g<sup>00</sup>)</code>,
whose covariant form is <code>n<sub>&mu;</sub> = (&minus;&alpha;, 0, 0, 0)</code> with lapse
<code>&alpha; = 1/&radic;(&minus;g<sup>00</sup>)</code>. In a stationary axisymmetric chart this is the
zero-angular-momentum observer, which stays timelike inside an ergosphere where <code>&part;<sub>t</sub></code>
is spacelike; it fails only where <code>g<sup>00</sup> &ge; 0</code> (on or inside a horizon of a
Boyer-Lindquist-type chart). For Schwarzschild it is the static observer and the spatial legs are the
radial, polar and azimuthal unit vectors. See
<a href=\"modelica://GeneralRelativity.Frames.observerTetrad\">observerTetrad</a> for the Gram-Schmidt step.</p>
</html>"));
end orthonormalTetrad;
