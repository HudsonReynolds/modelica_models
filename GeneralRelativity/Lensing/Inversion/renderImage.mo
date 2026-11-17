within GeneralRelativity.Lensing.Inversion;
function renderImage "Simulate the lensed image of Sersic sources on a pixel grid: the thin-lens equation per sub-pixel, then an optional Gaussian point-spread function"
  input Real u[:](each unit = "rad") "First image coordinate of the pixel centres (uniform spacing)";
  input Real v[:](each unit = "rad") "Second image coordinate of the pixel centres (uniform spacing)";
  input Real lens[7] "Lens parameters";
  input Real sources[:, 7] "One Sersic source per row";
  input Boolean pointMass = false;
  input Integer subpixels = 2 "Sub-pixel samples per axis (2 gives 4 per pixel)";
  input Real psfSigma(unit = "rad") = 0.0 "Gaussian point-spread function width; 0 for none";
  output Real image[size(u, 1), size(v, 1)] "Surface brightness at (u[i], v[j])";
protected
  Integer nu = size(u, 1);
  Integer nv = size(v, 1);
  Integer ns = size(sources, 1);
  Integer sub = max(subpixels, 1);
  Real du = if nu > 1 then u[2] - u[1] else 0.0;
  Real dv = if nv > 1 then v[2] - v[1] else 0.0;
  Integer half = if psfSigma > 0.0 and nu > 1 then integer(ceil(3.0 * psfSigma / abs(u[2] - u[1]))) else 0;
  Real weight = 1.0 / sub^2;
  // The lens, unpacked once (deflectionAt does the same per call).
  Real b = abs(lens[1]);
  Real lq = min(max(lens[4], 1.0e-3), 1.0);
  Real lc = cos(lens[5]);
  Real ls = sin(lens[5]);
  Real lsq = sqrt(max(1.0 - lq^2, 0.0));
  Real g1 = lens[6] * cos(2.0 * lens[7]);
  Real g2 = lens[6] * sin(2.0 * lens[7]);
  // The sources, unpacked once (sersicSource does the same per call): position, brightness, 1/re, axis frame, 1/q, 1/n, b_n.
  Real sx[size(sources, 1)];
  Real sy[size(sources, 1)];
  Real sIe[size(sources, 1)];
  Real sInvRe[size(sources, 1)];
  Real sc[size(sources, 1)];
  Real ss[size(sources, 1)];
  Real sInvQ[size(sources, 1)];
  Real sInvN[size(sources, 1)];
  Real sBn[size(sources, 1)];
  Real n;
  Real x;
  Real y;
  Real dx;
  Real dy;
  Real wx;
  Real wy;
  Real psi;
  Real ax;
  Real ay;
  Real awx;
  Real awy;
  Real r2;
  Real bx;
  Real by;
  Real xp;
  Real yp;
  Real R;
  Real acc;
  Real kernel[2 * (if psfSigma > 0.0 and size(u, 1) > 1 then integer(ceil(3.0 * psfSigma / abs(u[2] - u[1]))) else 0) + 1];
  Real blurred[size(u, 1), size(v, 1)];
  Real norm;
  Integer k;
algorithm
  for s in 1:ns loop
    sx[s] := sources[s, 1];
    sy[s] := sources[s, 2];
    sIe[s] := sources[s, 3];
    sInvRe[s] := 1.0 / max(sources[s, 4], 1.0e-30);
    // An axis ratio above 1 is the same ellipse with the axes swapped (sersicSource).
    if sources[s, 5] > 1.0 then
      sInvQ[s] := sources[s, 5];
      sc[s] := -sin(sources[s, 6]);
      ss[s] := cos(sources[s, 6]);
    else
      sInvQ[s] := 1.0 / max(sources[s, 5], 1.0e-3);
      sc[s] := cos(sources[s, 6]);
      ss[s] := sin(sources[s, 6]);
    end if;
    n := max(sources[s, 7], 0.2);
    sInvN[s] := 1.0 / n;
    sBn[s] := 2.0 * n - 1.0 / 3.0 + 4.0 / (405.0 * n) + 46.0 / (25515.0 * n^2);
  end for;
  for i in 1:nu loop
    for j in 1:nv loop
      acc := 0.0;
      for a in 1:sub loop
        for c in 1:sub loop
          x := u[i] + ((a - 0.5) / sub - 0.5) * du;
          y := v[j] + ((c - 0.5) / sub - 0.5) * dv;
          dx := x - lens[2];
          dy := y - lens[3];
          // Deflection: point mass or singular isothermal ellipse (major axis along the lens frame's x), then external shear.
          if pointMass then
            r2 := dx * dx + dy * dy;
            if r2 > 0.0 then
              ax := b * b * dx / r2;
              ay := b * b * dy / r2;
            else
              ax := 0.0;
              ay := 0.0;
            end if;
          else
            wx := lc * dx + ls * dy;
            wy := -ls * dx + lc * dy;
            psi := sqrt(lq * lq * wx * wx + wy * wy);
            if psi <= 0.0 then
              awx := 0.0;
              awy := 0.0;
            elseif lsq < 1.0e-6 then
              awx := b * wx / psi;
              awy := b * wy / psi;
            else
              awx := b * lq / lsq * atan(lsq * wx / psi);
              awy := b * lq / lsq * 0.5 * log((psi + lsq * wy) / (psi - lsq * wy));
            end if;
            ax := lc * awx - ls * awy;
            ay := ls * awx + lc * awy;
          end if;
          bx := x - ax - (g1 * dx + g2 * dy);
          by := y - ay - (g2 * dx - g1 * dy);
          for s in 1:ns loop
            xp := sc[s] * (bx - sx[s]) + ss[s] * (by - sy[s]);
            yp := (-ss[s] * (bx - sx[s]) + sc[s] * (by - sy[s])) * sInvQ[s];
            R := sqrt(xp * xp + yp * yp) * sInvRe[s];
            acc := acc + sIe[s] * exp(-sBn[s] * ((if sInvN[s] == 1.0 then R else R^sInvN[s]) - 1.0));
          end for;
        end for;
      end for;
      image[i, j] := weight * acc;
    end for;
  end for;
  // Separable Gaussian blur, kernel truncated at 3 sigma, renormalized at the edges.
  if half > 0 then
    for m in -half:half loop
      kernel[m + half + 1] := exp(-0.5 * (m * du / psfSigma)^2);
    end for;
    for i in 1:nu loop
      for j in 1:nv loop
        acc := 0.0;
        norm := 0.0;
        for m in -half:half loop
          k := i + m;
          if k >= 1 and k <= nu then
            acc := acc + kernel[m + half + 1] * image[k, j];
            norm := norm + kernel[m + half + 1];
          end if;
        end for;
        blurred[i, j] := acc / norm;
      end for;
    end for;
    for i in 1:nu loop
      for j in 1:nv loop
        acc := 0.0;
        norm := 0.0;
        for m in -half:half loop
          k := j + m;
          if k >= 1 and k <= nv then
            acc := acc + kernel[m + half + 1] * blurred[i, k];
            norm := norm + kernel[m + half + 1];
          end if;
        end for;
        image[i, j] := acc / norm;
      end for;
    end for;
  end if;
  annotation(Documentation(info="<html>
<p>Each sub-pixel at &theta; is traced to the source plane, <code>&beta; = &theta; &minus; &alpha;(&theta;)</code>, and shows the
summed source brightness there; arcs, rings and multiple images appear by themselves because lensing conserves surface
brightness. The point-spread function uses the spacing of <code>u</code> as the pixel size. This is the same calculation as
the viewer's lensing shader, on a grid small enough to fit repeatedly.</p>
<p>The loop is written with scalars only, the lens and source constants unpacked once, because the fit renders the image
dozens of times per iteration. It computes exactly what
<a href=\"modelica://GeneralRelativity.Lensing.Inversion.deflectionAt\">deflectionAt</a> and
<a href=\"modelica://GeneralRelativity.Lensing.Inversion.sersicSource\">sersicSource</a> compute, and the tests check the two
against each other.</p>
</html>"));
end renderImage;
