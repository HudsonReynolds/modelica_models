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
  Real theta[2];
  Real beta[2];
  Real acc;
  Real weight = 1.0 / sub^2;
  Real kernel[2 * (if psfSigma > 0.0 and size(u, 1) > 1 then integer(ceil(3.0 * psfSigma / abs(u[2] - u[1]))) else 0) + 1];
  Real blurred[size(u, 1), size(v, 1)];
  Real norm;
  Integer k;
algorithm
  for i in 1:nu loop
    for j in 1:nv loop
      acc := 0.0;
      for a in 1:sub loop
        for b in 1:sub loop
          theta := {u[i] + ((a - 0.5) / sub - 0.5) * du, v[j] + ((b - 0.5) / sub - 0.5) * dv};
          beta := theta - deflectionAt(theta, lens, pointMass);
          for s in 1:ns loop
            acc := acc + sersicSource(beta, sources[s, :]);
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
<p>Each sub-pixel at &theta; is traced to the source plane, <code>&beta; = &theta; &minus; &alpha;(&theta;)</code>
(<a href=\"modelica://GeneralRelativity.Lensing.Inversion.deflectionAt\">deflectionAt</a>), and shows the summed source
brightness there (<a href=\"modelica://GeneralRelativity.Lensing.Inversion.sersicSource\">sersicSource</a>); arcs, rings and
multiple images appear by themselves because lensing conserves surface brightness. The point-spread function uses the
spacing of <code>u</code> as the pixel size. This is the same calculation as the viewer's lensing shader, on a grid small enough
to fit repeatedly.</p>
</html>"));
end renderImage;
