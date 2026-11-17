within GeneralRelativity.Lensing.Inversion;
function sersicSource "Surface brightness of an elliptical Sersic source at a source-plane position"
  input Real beta[2](each unit = "rad") "Source-plane position";
  input Real source[7] "x, y, Ie, re, q, position angle, Sersic index n (Inversion package documentation)";
  output Real I "Surface brightness";
protected
  // An axis ratio above 1 is the same ellipse with the axes swapped: q -> 1/q, position angle + 90 degrees.
  // Letting q cross 1 keeps the fit smooth instead of pinning it at a bound.
  Boolean swapped = source[5] > 1.0;
  Real q = if swapped then 1.0 / source[5] else max(source[5], 1.0e-3);
  Real c = if swapped then -sin(source[6]) else cos(source[6]);
  Real s = if swapped then cos(source[6]) else sin(source[6]);
  Real dx = beta[1] - source[1];
  Real dy = beta[2] - source[2];
  Real xp = c * dx + s * dy "Along the major axis";
  Real yp = (-s * dx + c * dy) / q "Along the minor axis, stretched to a circle";
  Real n = max(source[7], 0.2);
  Real bn = 2.0 * n - 1.0 / 3.0 + 4.0 / (405.0 * n) + 46.0 / (25515.0 * n^2) "Ciotti and Bertin (1999)";
  Real R = sqrt(xp^2 + yp^2) / max(source[4], 1.0e-30);
algorithm
  I := source[3] * exp(-bn * (R^(1.0 / n) - 1.0));
  annotation(Inline = true, Documentation(info="<html>
<p><code>I(R) = I<sub>e</sub> exp(&minus;b<sub>n</sub>[(R/r<sub>e</sub>)<sup>1/n</sup> &minus; 1])</code> with the elliptical radius
<code>R = &radic;(x'&sup2; + y'&sup2;/q&sup2;)</code> in the frame of the major axis; <code>b<sub>n</sub></code> makes <code>r<sub>e</sub></code>
the half-light radius. Lensing conserves surface brightness, so an image pixel at &theta; simply shows this at
<code>&beta; = &theta; &minus; &alpha;(&theta;)</code>.</p>
</html>"));
end sersicSource;
