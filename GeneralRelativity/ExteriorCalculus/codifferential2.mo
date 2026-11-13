within GeneralRelativity.ExteriorCalculus;
function codifferential2 "Codifferential of a 2-form field, (delta F)_b = -nabla^a F_ab = (*d*F)_b"
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric function";
  input TwoFormField field "2-form field";
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Metric parameters";
  input Real q[:] "Field parameters";
  input Real step[4] "Finite-difference step per coordinate";
  output Real w[4] "1-form (delta F)_b = -g_bd (1/sqrt|g|) d_a (sqrt|g| F^ad)";
protected
  Real xShift[4];
  Real gShift[4, 4];
  Real gInv[4, 4];
  Boolean accepted;
  Real FShift[4, 4];
  Real densityPlus[4, 4];
  Real densityMinus[4, 4];
  Real flux[4] "d_a (sqrt|g| F^ad)";
  Real g[4, 4];
algorithm
  flux := zeros(4);
  for a in 1:4 loop
    xShift := x;
    xShift[a] := x[a] + step[a];
    gShift := metric(xShift, p);
    FShift := field(xShift, q);
    (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(gShift);
    densityPlus := volumeElement(gShift) * (gInv * FShift * gInv);
    xShift[a] := x[a] - step[a];
    gShift := metric(xShift, p);
    FShift := field(xShift, q);
    (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(gShift);
    densityMinus := volumeElement(gShift) * (gInv * FShift * gInv);
    flux := flux + (densityPlus[a, :] - densityMinus[a, :]) / (2.0 * step[a]);
  end for;
  g := metric(x, p);
  w := -(g * flux) / volumeElement(g);
  annotation(Documentation(info="<html>
<p>Uses <code>&nabla;<sub>a</sub>F<sup>ab</sup> = (1/&radic;|g|)&part;<sub>a</sub>(&radic;|g|F<sup>ab</sup>)</code>
for antisymmetric <code>F</code>. The inhomogeneous Maxwell equations read
<code>&nabla;<sub>a</sub>F<sup>ab</sup> = -&mu;<sub>0</sub>J<sup>b</sup></code> or
<code>+&mu;<sub>0</sub>J<sup>b</sup></code> depending on the sign of <code>F</code>;
in vacuum <code>&delta;F = 0</code> either way. See
<a href=\"modelica://GeneralRelativity.ExteriorCalculus.codifferential1\">codifferential1</a>
for the sign of <code>&delta;</code>.</p>
</html>"));
end codifferential2;
