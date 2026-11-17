within GeneralRelativity.Lensing.Inversion;
function invert "Steps 2 to 5 of a lens inversion: read the ring off the image, start short fits from two candidate lens positions, and refine the better one until the simulation matches"
  input Real u[:](each unit = "rad");
  input Real v[:](each unit = "rad");
  input Real observed[size(u, 1), size(v, 1)];
  input Real noiseSigma "Pixel noise standard deviation";
  input Boolean pointMass = false;
  input Real psfSigma(unit = "rad") = 0.0;
  input Integer maxIterations = 30 "Iterations of the final refinement";
  input Integer startIterations = 4 "Iterations of each trial start";
  input Real threshold = 0.2 "Fraction of the peak below which pixels do not count for the ring";
  input Integer subpixels = 2;
  output Real ringRadius(unit = "rad") "Radius of the circle through the arcs: the radius-based Einstein radius";
  output Real ringCentre[2](each unit = "rad");
  output Real lens0[7] "The start that won, before any iteration";
  output Real sources0[1, 7];
  output Real lens[7] "Fitted lens";
  output Real sources[1, 7] "Fitted source";
  output Real chi2History[startIterations + maxIterations + 1] "Cost from the winning start through the refinement, padded with the final value";
  output Real lensHistory[startIterations + maxIterations + 1, 7];
  output Integer iterations "Iterations of the winning start plus the refinement";
  output Boolean converged;
  output Real lensSigma[7];
  output Real sourcesSigma[1, 7];
  output Real lensCovariance[7, 7];
protected
  Real pixel = if size(u, 1) > 1 then abs(u[2] - u[1]) else 1.0e-9;
  Real spread;
  Real light[2] "Light centroid of the arcs";
  Real starts[8, 2] "Candidate lens positions";
  Integer nStarts;
  Real shapes[3, 2] "Candidate (axis ratio, position angle): nearly round, and from the second moments of the arcs either way round";
  Real sourceRadius0 "Starting source size: the radial thickness of the arcs, less the seeing";
  Real peak = max(observed);
  Real w;
  Real dx;
  Real dy;
  Real sumW;
  Real mxx;
  Real myy;
  Real mxy;
  Real lambdaMax;
  Real lambdaMin;
  Real trialCost;
  Real trialLens[7];
  Real trialSources[1, 7];
  Real trialChi2[startIterations + 1];
  Real trialHistory[startIterations + 1, 7];
  Integer trialIterations;
  Boolean trialConverged;
  Real trialLensSigma[7];
  Real trialSourcesSigma[1, 7];
  Real trialCovariance[7, 7];
  Real bestChi2 = 0.0;
  Real bestLens[7];
  Real bestSources[1, 7];
  Real bestChi2History[startIterations + 1];
  Real bestLensHistory[startIterations + 1, 7];
  Integer bestIterations = 0;
  Real lensStart[7];
  Real sourcesStart[1, 7];
  Real refineChi2[maxIterations + 1];
  Real refineHistory[maxIterations + 1, 7];
  Integer refineIterations;
  Integer k;
algorithm
  // Step 2: the circle through the bright pixels. Its radius is the equal-area Einstein radius (einsteinMass).
  (ringRadius, ringCentre, spread) := einsteinRadiusFromImage(u, v, observed, threshold);
  light := brightnessCentroid(u, v, observed, threshold);
  // For an isothermal lens the images of a source at beta lie on a circle of radius thetaE about beta, and their light centroid
  // is at 2 beta, so the lens is near 2 * circle centre - light centroid. For a ring both coincide with the lens; when the images
  // are not those of a round lens, the circle centre is the better guess. Both are tried.
  starts[1, :] := ringCentre;
  starts[2, :] := 2.0 * ringCentre - light;
  // Six more around the circle centre, at a third of the ring radius: where off-axis sources put the lens otherwise.
  for m in 1:6 loop
    starts[2 + m, :] := ringCentre + 0.35 * ringRadius * {cos(GeneralRelativity.Constants.pi * (m - 1) / 3.0), sin(GeneralRelativity.Constants.pi * (m - 1) / 3.0)};
  end for;
  nStarts := 8;
  if sqrt((starts[2, 1] - starts[1, 1])^2 + (starts[2, 2] - starts[1, 2])^2) <= pixel then
    starts[2, :] := starts[8, :];
    nStarts := 7;
  end if;
  // An isothermal lens hardly stretches images radially, so the arcs are about as thick as the source convolved with the seeing.
  // The half-light radius of a profile with standard deviation s is about 1.2 s.
  sourceRadius0 := min(max(1.2 * sqrt(max(spread^2 - psfSigma^2, 0.0)), pixel), 0.3 * ringRadius);
  // The source starts at the circle centre: where an isothermal lens puts it.
  sourcesStart := {{ringCentre[1], ringCentre[2], peak, sourceRadius0, 0.9, 0.0, 1.0}};
  bestIterations := 0;
  for s in 1:nStarts loop
    // The arcs trace the critical curve, an ellipse with semi-axes b and q b along the lens's major axis: the second moments of
    // the bright pixels about the candidate centre give its orientation and, as the square root of their ratio, q.
    sumW := 0.0;
    mxx := 0.0;
    myy := 0.0;
    mxy := 0.0;
    for i in 1:size(u, 1) loop
      for j in 1:size(v, 1) loop
        w := observed[i, j] - threshold * peak;
        if w > 0.0 then
          dx := u[i] - starts[s, 1];
          dy := v[j] - starts[s, 2];
          sumW := sumW + w;
          mxx := mxx + w * dx * dx;
          myy := myy + w * dy * dy;
          mxy := mxy + w * dx * dy;
        end if;
      end for;
    end for;
    lambdaMax := 0.5 * (mxx + myy) + sqrt(0.25 * (mxx - myy)^2 + mxy^2);
    lambdaMin := 0.5 * (mxx + myy) - sqrt(0.25 * (mxx - myy)^2 + mxy^2);
    shapes[1, :] := {0.9, 0.0};
    shapes[2, :] := {min(max(sqrt(max(lambdaMin, 0.0) / max(lambdaMax, 1.0e-300)), 0.4), 0.95), 0.5 * atan2(2.0 * mxy, mxx - myy)};
    // Which images are brightest decides whether the moments line up with the major or the minor axis, so try both.
    shapes[3, :] := {shapes[2, 1], shapes[2, 2] + 0.5 * GeneralRelativity.Constants.pi};
    for shape in 1:(if pointMass then 1 else 3) loop
      lensStart := {if pointMass then ringRadius else ringRadius / sqrt(shapes[shape, 1]), starts[s, 1], starts[s, 2], shapes[shape, 1], shapes[shape, 2], 0.0, 0.0};
      (trialLens, trialSources, trialChi2, trialHistory, trialIterations, trialConverged, trialLensSigma, trialSourcesSigma, trialCovariance) :=
        fitLens(u, v, observed, noiseSigma, lensStart, sourcesStart, pointMass = pointMass, subpixels = subpixels, psfSigma = psfSigma, maxIterations = startIterations);
      trialCost := trialChi2[startIterations + 1];
      if (s == 1 and shape == 1) or trialCost < bestChi2 then
        bestChi2 := trialCost;
        bestLens := trialLens;
        bestSources := trialSources;
        bestChi2History := trialChi2;
        bestLensHistory := trialHistory;
        bestIterations := trialIterations;
        lens0 := lensStart;
        sources0 := sourcesStart;
      end if;
    end for;
  end for;
  // Step 5: refine the better start to convergence.
  (lens, sources, refineChi2, refineHistory, refineIterations, converged, lensSigma, sourcesSigma, lensCovariance) :=
    fitLens(u, v, observed, noiseSigma, bestLens, bestSources, pointMass = pointMass, subpixels = subpixels, psfSigma = psfSigma, maxIterations = maxIterations);
  iterations := bestIterations + refineIterations;
  k := 0;
  for i in 1:bestIterations + 1 loop
    k := k + 1;
    chi2History[k] := bestChi2History[i];
    lensHistory[k, :] := bestLensHistory[i, :];
  end for;
  for i in 2:maxIterations + 1 loop
    k := k + 1;
    chi2History[k] := refineChi2[i];
    lensHistory[k, :] := refineHistory[i, :];
  end for;
  for i in k + 1:startIterations + maxIterations + 1 loop
    chi2History[i] := refineChi2[maxIterations + 1];
    lensHistory[i, :] := refineHistory[maxIterations + 1, :];
  end for;
  annotation(Documentation(info="<html>
<p>The usual way to call the inversion. Pixel fits of lenses have local minima, so it starts from several measured guesses:</p>
<ul>
<li><b>Position.</b> A fit started at the centre of the circle through the arcs fails for an off-axis source: that circle is
centred on the <i>source</i>, and a lens placed there with the circle's radius explains the image positions but not their
brightness. Brightness breaks the tie: for an isothermal lens the light centroid of the images sits twice as far from the lens as
the circle centre does. The circle centre, <code>2 &times; centre &minus; centroid</code> and six points a third of the ring radius
around the centre are tried.</li>
<li><b>Shape.</b> The arcs trace the critical curve, an ellipse with semi-axes <code>b</code> and <code>qb</code>, so the second
moments of the bright pixels give the position angle and q. Which images are brightest decides whether the moments line up with
the major or the minor axis, so a nearly round start (q = 0.9) and the moment-based shape at both orientations are tried; a point mass
has no shape.</li>
<li><b>Source size.</b> An isothermal lens hardly magnifies radially, so the radial thickness of the arcs (less the seeing) is the
source size. Starting with a large source instead converges to a smooth ring when the truth is a compact quasar-like source.</li>
</ul>
<p><code>invert</code> runs <a href=\"modelica://GeneralRelativity.Lensing.Inversion.fitLens\">fitLens</a> for
<code>startIterations</code> from each start, keeps the lowest &chi;&sup2;, and refines it for up to <code>maxIterations</code> more. The histories run through
the winning start and the refinement, so their first entry is the cost of the first estimate.</p>
<p>Every start has the circle's radius as its equal-area Einstein radius (<code>b = R/&radic;q<sub>0</sub></code>), no shear, and an
exponential source at the circle centre with the peak brightness. Use
<a href=\"modelica://GeneralRelativity.Lensing.Inversion.fitLens\">fitLens</a> directly to supply your own start or priors.</p>
</html>"));
end invert;
