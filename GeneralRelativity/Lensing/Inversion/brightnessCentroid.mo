within GeneralRelativity.Lensing.Inversion;
function brightnessCentroid "Flux-weighted centroid of the pixels above a fraction of the peak: the lens position for a ring or a symmetric pair of arcs"
  input Real u[:](each unit = "rad");
  input Real v[:](each unit = "rad");
  input Real image[size(u, 1), size(v, 1)];
  input Real threshold = 0.2 "Fraction of the peak brightness below which pixels are ignored (noise)";
  output Real centre[2](each unit = "rad");
protected
  Real peak = max(image);
  Real w;
  Real sumW = 0.0;
  Real sumX = 0.0;
  Real sumY = 0.0;
algorithm
  for i in 1:size(u, 1) loop
    for j in 1:size(v, 1) loop
      w := image[i, j] - threshold * peak;
      if w > 0.0 then
        sumW := sumW + w;
        sumX := sumX + w * u[i];
        sumY := sumY + w * v[j];
      end if;
    end for;
  end for;
  centre := if sumW > 0.0 then {sumX, sumY} / sumW else {0.0, 0.0};
end brightnessCentroid;
