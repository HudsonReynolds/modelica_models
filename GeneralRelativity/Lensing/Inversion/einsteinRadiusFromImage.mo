within GeneralRelativity.Lensing.Inversion;
function einsteinRadiusFromImage "Einstein radius read off an image: the circle fitted through the bright pixels, weighted by brightness"
  input Real u[:](each unit = "rad");
  input Real v[:](each unit = "rad");
  input Real image[size(u, 1), size(v, 1)];
  input Real threshold = 0.2 "Fraction of the peak brightness below which pixels are ignored (noise)";
  output Real thetaE(unit = "rad") "Radius of the fitted circle: the Einstein radius";
  output Real centre[2](each unit = "rad") "Centre of the fitted circle: the lens position, shifted by the source offset when the arcs are one-sided";
  output Real spread(unit = "rad") "Flux-weighted standard deviation of the pixel distances from the circle: the ring's width, or how badly a circle fits";
protected
  Real peak = max(image);
  Real scale = max(max(abs(u)), max(abs(v))) "Coordinates are fitted in units of the field so the normal equations are well scaled";
  Real w;
  Real x;
  Real y;
  Real r;
  Real A[3, 3] = zeros(3, 3);
  Real b[3, 1] = zeros(3, 1);
  Real coefficients[3, 1];
  Boolean ok;
  Real sumW = 0.0;
  Real sumD = 0.0;
  Real sumD2 = 0.0;
algorithm
  // Algebraic (Kasa) circle fit: minimize sum w (x^2 + y^2 + D x + E y + F)^2, linear in D, E, F.
  for i in 1:size(u, 1) loop
    for j in 1:size(v, 1) loop
      w := image[i, j] - threshold * peak;
      if w > 0.0 then
        x := u[i] / scale;
        y := v[j] / scale;
        A := A + w * [x * x, x * y, x; x * y, y * y, y; x, y, 1.0];
        b := b - w * (x * x + y * y) * [x; y; 1.0];
      end if;
    end for;
  end for;
  (coefficients, ok) := LinearAlgebra.solve(A, b);
  if ok then
    centre := {-0.5 * coefficients[1, 1], -0.5 * coefficients[2, 1]};
    thetaE := sqrt(max(centre[1]^2 + centre[2]^2 - coefficients[3, 1], 0.0)) * scale;
    centre := centre * scale;
  else
    // Degenerate (too few bright pixels): fall back to the centroid and the mean radius about it.
    centre := brightnessCentroid(u, v, image, threshold);
    thetaE := 0.0;
    for i in 1:size(u, 1) loop
      for j in 1:size(v, 1) loop
        w := image[i, j] - threshold * peak;
        if w > 0.0 then
          sumW := sumW + w;
          thetaE := thetaE + w * sqrt((u[i] - centre[1])^2 + (v[j] - centre[2])^2);
        end if;
      end for;
    end for;
    thetaE := if sumW > 0.0 then thetaE / sumW else 0.0;
    sumW := 0.0;
  end if;
  for i in 1:size(u, 1) loop
    for j in 1:size(v, 1) loop
      w := image[i, j] - threshold * peak;
      if w > 0.0 then
        r := sqrt((u[i] - centre[1])^2 + (v[j] - centre[2])^2) - thetaE;
        sumW := sumW + w;
        sumD := sumD + w * r;
        sumD2 := sumD2 + w * r^2;
      end if;
    end for;
  end for;
  spread := if sumW > 0.0 then sqrt(max(sumD2 / sumW - (sumD / sumW)^2, 0.0)) else 0.0;
  annotation(Documentation(info="<html>
<p>A source on the axis of a symmetric lens makes a ring of radius &theta;<sub>E</sub>. Off the axis an isothermal sphere
makes two images at <code>&theta;<sub>E</sub> &plusmn; &beta;</code> on opposite sides of the lens, which lie on a circle of radius
&theta;<sub>E</sub> centred at the lens plus the source offset, so a circle fitted to the arcs recovers &theta;<sub>E</sub> exactly for
that lens and approximately for an ellipse, a point mass or an extended source. The fit is the algebraic (K&aring;sa) least-squares
circle through the pixels above the threshold, each weighted by its brightness above it; with fewer than three independent bright
pixels it falls back to the brightness centroid and the mean radius. It is the starting point for
<a href=\"modelica://GeneralRelativity.Lensing.Inversion.fitLens\">fitLens</a>, not the answer; the mass it implies is
<a href=\"modelica://GeneralRelativity.Lensing.massInEinsteinRadius\">massInEinsteinRadius</a>.</p>
</html>"));
end einsteinRadiusFromImage;
