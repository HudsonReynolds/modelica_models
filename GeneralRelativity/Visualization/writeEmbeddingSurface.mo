within GeneralRelativity.Visualization;
impure function writeEmbeddingSurface "Append the equatorial embedding surface of a stationary axisymmetric metric, coloured by the Kretschmann scalar"
  input String fileName;
  input GeneralRelativity.Interfaces.MetricFunction metric "Metric in {ct, r, theta, phi}";
  input Real p[:] "Metric parameters";
  input Real step[4] "Finite-difference step for curvature";
  input Real rMin "Innermost radius, outside any horizon";
  input Real rMax "Outermost radius";
  input Integer nR = 40 "Radial rows";
  input Integer nPhi = 64 "Azimuthal samples per row";
protected
  Real r;
  Real rPrevious;
  Real slope;
  Real slopePrevious;
  Real z;
  Real x[4];
  Real g[4, 4];
  Real gInv[4, 4];
  Boolean accepted;
  Real dg[4, 4, 4];
  Real ddg[4, 4, 4, 4];
  Real Gamma[4, 4, 4];
  Real riemann[4, 4, 4, 4];
  Real kretschmann;
  Real phi;
  Real circumferentialRadius "R = sqrt(g_phiphi), the proper circumference over 2 pi";
  String row;
algorithm
  z := 0.0;
  rPrevious := rMin;
  slopePrevious := 0.0;
  for i in 1:nR loop
    r := rMin + (rMax - rMin) * (i - 1) / (nR - 1);
    x := {0.0, r, GeneralRelativity.Constants.pi / 2.0, 0.0};
    g := metric(x, p);
    (gInv, accepted) := GeneralRelativity.Tensors.inverseMetric(g);
    dg := GeneralRelativity.Metrics.derivative(metric, x, p, step);
    ddg := GeneralRelativity.Metrics.secondDerivative(metric, x, p, step);
    Gamma := GeneralRelativity.Tensors.christoffel(gInv, dg);
    riemann := GeneralRelativity.Tensors.riemannLower(g, ddg, Gamma);
    kretschmann := GeneralRelativity.Tensors.kretschmann(gInv, riemann);
    // Equatorial geometry ds^2 = g_rr dr^2 + g_phiphi dphi^2 embedded as dz^2 + dR^2 + R^2 dphi^2:
    // R = sqrt(g_phiphi) and dz/dr = sqrt(g_rr - (dR/dr)^2), integrated by the trapezoidal rule from rMin.
    // For Schwarzschild R = r and this is Flamm's sqrt(g_rr - 1).
    circumferentialRadius := sqrt(g[4, 4]);
    slope := sqrt(max(g[2, 2] - (dg[2, 4, 4] / (2.0 * circumferentialRadius))^2, 0.0));
    if i > 1 then
      z := z + 0.5 * (slope + slopePrevious) * (r - rPrevious);
    end if;
    row := "{\"type\":\"surfaceRow\",\"r\":" + jsonNumber(r) + ",\"R\":" + jsonNumber(circumferentialRadius) + ",\"z\":" + jsonNumber(z)
      + ",\"kretschmann\":" + jsonNumber(kretschmann) + ",\"nPhi\":" + String(nPhi) + "}";
    Modelica.Utilities.Streams.print(row, fileName);
    rPrevious := r;
    slopePrevious := slope;
  end for;
  annotation(Documentation(info="<html>
<p>The equatorial plane <code>t = const, &theta; = &pi;/2</code> is embedded in Euclidean space as a
surface of revolution with circumferential radius <code>R = &radic;g<sub>&phi;&phi;</sub></code> and height
<code>dz/dr = &radic;(g<sub>rr</sub> &minus; (dR/dr)&sup2;)</code>, exact for any stationary axisymmetric metric
whose equatorial plane has no <code>g<sub>r&phi;</sub></code> term (Schwarzschild, Reissner-Nordstr&ouml;m,
Kerr); for Schwarzschild it is Flamm's paraboloid. Where <code>(dR/dr)&sup2; &gt; g<sub>rr</sub></code> the
plane does not embed and the slope is clipped to zero.</p>
<p>Each row is one radius; the viewer revolves it into <code>nPhi</code> vertices at radius <code>R</code>.
The height starts at zero at <code>rMin</code>. The trapezoidal integration of
<code>sqrt(g<sub>rr</sub> - 1)</code> is second order except near a horizon where
the slope diverges like <code>(r - r<sub>s</sub>)<sup>-1/2</sup></code>; start a row or
two outside it.</p>
</html>"));
end writeEmbeddingSurface;
