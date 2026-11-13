within GeneralRelativity.ExteriorCalculus;
function codifferential1 "Codifferential of a 1-form field, delta w = -nabla^a w_a = *d*w"
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric function";
  input OneFormField field "1-form field";
  input Real x[4] "Coordinates; x[1] = c*t";
  input Real p[:] "Metric parameters";
  input Real q[:] "Field parameters";
  input Real step[4] "Finite-difference step per coordinate";
  output Real f "0-form delta w = -(1/sqrt|g|) d_a (sqrt|g| g^ab w_b)";
protected
  Real xShift[4];
  Real gShift[4, 4];
  Real gInv[4, 4];
  Boolean accepted;
  Real wShift[4];
  Real densityPlus;
  Real densityMinus;
  Real flux;
algorithm
  flux := 0.0;
  for a in 1:4 loop
    xShift := x;
    xShift[a] := x[a] + step[a];
    gShift := metric(xShift, p);
    wShift := field(xShift, q);
    (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(gShift);
    densityPlus := volumeElement(gShift) * (gInv[a, :] * wShift);
    xShift[a] := x[a] - step[a];
    gShift := metric(xShift, p);
    wShift := field(xShift, q);
    (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(gShift);
    densityMinus := volumeElement(gShift) * (gInv[a, :] * wShift);
    flux := flux + (densityPlus - densityMinus) / (2.0 * step[a]);
  end for;
  gShift := metric(x, p);
  f := -flux / volumeElement(gShift);
  annotation(Documentation(info="<html>
<p>Sign convention: <code>&delta; = -div</code>, the formal adjoint of <code>d</code>
for the indefinite inner product. In four Lorentzian dimensions with the Hodge
dual of this package, <code>&delta; = *d*</code> on forms of every degree (the
general rule is <code>&delta; = (-1)<sup>n(p+1)+1</sup> s *d*</code> with
<code>s = -1</code> the sign of <code>det g</code>). The divergence is taken in
conservation form, so no Christoffel symbols are needed.</p>
</html>"));
end codifferential1;
