within GeneralRelativity.Perturbation;
function weakFieldMetric "Newtonian-limit metric g = eta + h with h_00 = h_ii = -2 Phi/c^2 and optional gravitomagnetic h_0i"
  input Real phi(unit = "m2/s2") "Newtonian potential, negative near mass";
  input Real h0i[3] = zeros(3) "Mixed time-space perturbation h_0i (dimensionless)";
  input Real c(unit = "m/s") = 1.0 "Speed of light";
  output Real g[4, 4] "Covariant metric";
protected
  Real h = -2.0 * phi / c^2;
algorithm
  g := zeros(4, 4);
  g[1, 1] := -1.0 + h;
  for i in 1:3 loop
    g[i + 1, i + 1] := 1.0 + h;
    g[1, i + 1] := h0i[i];
    g[i + 1, 1] := h0i[i];
  end for;
  annotation(Documentation(info="<html>
<p>Spatial coordinates are isotropic. With <code>&Phi; = -GM/r</code> this is
Schwarzschild in isotropic coordinates to first order. The redshift between two
static clocks is <code>&Delta;&Phi;/c&sup2;</code>; light bending and perihelion
precession are reproduced to leading order because the spatial part carries the
same <code>h</code> as the time part (<code>&gamma; = 1</code>).</p>
</html>"));
end weakFieldMetric;
