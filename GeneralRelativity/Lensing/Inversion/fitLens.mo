within GeneralRelativity.Lensing.Inversion;
function fitLens "Refine lens and sources until the simulated image matches the observed one: Levenberg-Marquardt least squares on the pixels, with optional Gaussian priors (the iterated extended Kalman update)"
  input Real u[:](each unit = "rad");
  input Real v[:](each unit = "rad");
  input Real observed[size(u, 1), size(v, 1)];
  input Real noiseSigma "Pixel noise standard deviation";
  input Real lens0[7] "Starting lens (Inversion package documentation)";
  input Real sources0[:, 7] "Starting sources, one per row";
  input Boolean freeLens[7] = fill(true, 7) "false fixes a lens parameter at its starting value";
  input Boolean freeSources[size(sources0, 1), 7] = fill(true, size(sources0, 1), 7);
  input Real priorSigmaLens[7] = zeros(7) "Prior standard deviation about lens0 per parameter; 0 means no prior";
  input Real priorSigmaSources[size(sources0, 1), 7] = zeros(size(sources0, 1), 7);
  input Boolean pointMass = false "true: point-mass lens (q and the position angle are then held fixed)";
  input Integer subpixels = 2;
  input Real psfSigma(unit = "rad") = 0.0;
  input Integer maxIterations = 30;
  input Real tolerance = 1.0e-6 "Stop when an accepted step lowers chi2 by less than this fraction";
  output Real lens[7] "Fitted lens";
  output Real sources[size(sources0, 1), 7] "Fitted sources";
  output Real chi2History[maxIterations + 1] "Cost before each iteration and after the last (padded with the final value)";
  output Real lensHistory[maxIterations + 1, 7] "Lens parameters alongside chi2History, for the mass estimate per iteration";
  output Integer iterations "Iterations taken";
  output Boolean converged "false only when the iteration limit stopped the fit";
  output Real lensSigma[7] "Standard error of each lens parameter from the curvature of chi2 at the fit; 0 when fixed";
  output Real sourcesSigma[size(sources0, 1), 7];
protected
  Integer ns = size(sources0, 1);
  Integer np = 7 + 7 * size(sources0, 1);
  Integer nu = size(u, 1);
  Integer nv = size(v, 1);
  Real pixel = if size(u, 1) > 1 then abs(u[2] - u[1]) else 1.0e-9;
  Real thetaScale = max(abs(lens0[1]), pixel);
  Real p[7 + 7 * size(sources0, 1)];
  Real p0[7 + 7 * size(sources0, 1)];
  Real pTry[7 + 7 * size(sources0, 1)];
  Real scale[7 + 7 * size(sources0, 1)];
  Real weight[7 + 7 * size(sources0, 1)] "Prior weight (scale / prior sigma)^2 in scaled coordinates";
  Boolean free[7 + 7 * size(sources0, 1)];
  Real simulated[size(u, 1), size(v, 1)];
  Real plus[size(u, 1), size(v, 1)];
  Real minus[size(u, 1), size(v, 1)];
  Real J[size(u, 1) * size(v, 1), 7 + 7 * size(sources0, 1)] "d residual / d scaled parameter";
  Real A[7 + 7 * size(sources0, 1), 7 + 7 * size(sources0, 1)];
  Real D[7 + 7 * size(sources0, 1), 7 + 7 * size(sources0, 1)];
  Real rhs[7 + 7 * size(sources0, 1), 1];
  Real delta[7 + 7 * size(sources0, 1), 1];
  Real C[7 + 7 * size(sources0, 1), 7 + 7 * size(sources0, 1)];
  Real srcTry[size(sources0, 1), 7];
  Real chi2;
  Real chi2Try;
  Real lambda = 1.0e-3;
  Real maxDiag;
  Real h = 1.0e-4 "Central-difference step in scaled units";
  Boolean accepted;
  Boolean ok;
  Boolean stop;
  Integer tries;
  Integer pix;
algorithm
  // Pack lens and sources into one vector, with the scale of each parameter, its freedom and its prior.
  p[1:7] := lens0;
  scale[1:7] := {thetaScale, thetaScale, thetaScale, 0.2, 0.3, 0.05, 0.3};
  for k in 1:7 loop
    free[k] := freeLens[k] and not (pointMass and (k == 4 or k == 5));
    weight[k] := if priorSigmaLens[k] > 0.0 then (scale[k] / priorSigmaLens[k])^2 else 0.0;
  end for;
  for s in 1:ns loop
    scale[7 + 7 * (s - 1) + 1:7 + 7 * (s - 1) + 7] := {thetaScale, thetaScale, max(abs(sources0[s, 3]), 1.0e-30), max(abs(sources0[s, 4]), 0.05 * thetaScale), 0.2, 0.3, 0.5};
    for k in 1:7 loop
      p[7 + 7 * (s - 1) + k] := sources0[s, k];
      free[7 + 7 * (s - 1) + k] := freeSources[s, k];
      weight[7 + 7 * (s - 1) + k] := if priorSigmaSources[s, k] > 0.0 then (scale[7 + 7 * (s - 1) + k] / priorSigmaSources[s, k])^2 else 0.0;
    end for;
  end for;
  p0 := p;
  // Cost at the start.
  srcTry := sources0;
  simulated := renderImage(u, v, p[1:7], srcTry, pointMass, subpixels, psfSigma);
  chi2 := 0.0;
  for i in 1:nu loop
    for j in 1:nv loop
      chi2 := chi2 + ((observed[i, j] - simulated[i, j]) / noiseSigma)^2;
    end for;
  end for;
  chi2History := fill(chi2, maxIterations + 1);
  for k in 1:maxIterations + 1 loop
    lensHistory[k, :] := p[1:7];
  end for;
  iterations := 0;
  converged := false;
  stop := maxIterations < 1;
  while not stop loop
    iterations := iterations + 1;
    // Jacobian of the weighted residual r = (observed - simulated) / sigma with respect to the scaled free parameters.
    J := zeros(nu * nv, np);
    for k in 1:np loop
      if free[k] then
        pTry := p;
        pTry[k] := p[k] + h * scale[k];
        for s in 1:ns loop
          srcTry[s, :] := pTry[7 + 7 * (s - 1) + 1:7 + 7 * (s - 1) + 7];
        end for;
        plus := renderImage(u, v, pTry[1:7], srcTry, pointMass, subpixels, psfSigma);
        pTry[k] := p[k] - h * scale[k];
        for s in 1:ns loop
          srcTry[s, :] := pTry[7 + 7 * (s - 1) + 1:7 + 7 * (s - 1) + 7];
        end for;
        minus := renderImage(u, v, pTry[1:7], srcTry, pointMass, subpixels, psfSigma);
        for i in 1:nu loop
          for j in 1:nv loop
            J[(i - 1) * nv + j, k] := -(plus[i, j] - minus[i, j]) / (2.0 * h * noiseSigma);
          end for;
        end for;
      end if;
    end for;
    // Normal equations (J'J + W + lambda D) delta = -(J' r + W (s - s0)); fixed parameters get a unit row and zero right-hand side.
    A := transpose(J) * J;
    rhs := zeros(np, 1);
    for i in 1:nu loop
      for j in 1:nv loop
        pix := (i - 1) * nv + j;
        for k in 1:np loop
          rhs[k, 1] := rhs[k, 1] - J[pix, k] * (observed[i, j] - simulated[i, j]) / noiseSigma;
        end for;
      end for;
    end for;
    maxDiag := 1.0e-300;
    for k in 1:np loop
      A[k, k] := A[k, k] + weight[k];
      rhs[k, 1] := rhs[k, 1] - weight[k] * (p[k] - p0[k]) / scale[k];
      if not free[k] then
        A[k, :] := zeros(np);
        A[:, k] := zeros(np);
        A[k, k] := 1.0;
        rhs[k, 1] := 0.0;
      end if;
      maxDiag := max(maxDiag, A[k, k]);
    end for;
    D := zeros(np, np);
    for k in 1:np loop
      D[k, k] := if free[k] then max(A[k, k], 1.0e-8 * maxDiag) else 0.0;
    end for;
    // Try damped steps until one lowers the cost.
    accepted := false;
    tries := 0;
    while not accepted and tries < 12 loop
      tries := tries + 1;
      (delta, ok) := LinearAlgebra.solve(A + lambda * D, rhs);
      if ok then
        for k in 1:np loop
          pTry[k] := p[k] + (if free[k] then delta[k, 1] * scale[k] else 0.0);
        end for;
        // Keep the parameters in their physical ranges.
        pTry[1] := abs(pTry[1]);
        pTry[4] := min(max(pTry[4], 0.1), 1.0);
        for s in 1:ns loop
          pTry[7 + 7 * (s - 1) + 3] := max(pTry[7 + 7 * (s - 1) + 3], 0.0);
          pTry[7 + 7 * (s - 1) + 4] := max(pTry[7 + 7 * (s - 1) + 4], 1.0e-3 * thetaScale);
          pTry[7 + 7 * (s - 1) + 5] := min(max(pTry[7 + 7 * (s - 1) + 5], 0.1), 10.0);
          pTry[7 + 7 * (s - 1) + 7] := min(max(pTry[7 + 7 * (s - 1) + 7], 0.3), 8.0);
          srcTry[s, :] := pTry[7 + 7 * (s - 1) + 1:7 + 7 * (s - 1) + 7];
        end for;
        plus := renderImage(u, v, pTry[1:7], srcTry, pointMass, subpixels, psfSigma);
        chi2Try := 0.0;
        for i in 1:nu loop
          for j in 1:nv loop
            chi2Try := chi2Try + ((observed[i, j] - plus[i, j]) / noiseSigma)^2;
          end for;
        end for;
        for k in 1:np loop
          chi2Try := chi2Try + weight[k] * ((pTry[k] - p0[k]) / scale[k])^2;
        end for;
        if chi2Try < chi2 then
          accepted := true;
          converged := chi2 - chi2Try < tolerance * chi2;
          p := pTry;
          simulated := plus;
          chi2 := chi2Try;
          lambda := max(lambda / 3.0, 1.0e-12);
        else
          lambda := lambda * 10.0;
        end if;
      else
        lambda := lambda * 10.0;
      end if;
    end while;
    chi2History[iterations + 1] := chi2;
    lensHistory[iterations + 1, :] := p[1:7];
    // A step that no damping can improve means the fit is at its numerical minimum.
    if not accepted then
      converged := true;
    end if;
    stop := converged or iterations >= maxIterations;
  end while;
  for k in iterations + 2:maxIterations + 1 loop
    chi2History[k] := chi2;
    lensHistory[k, :] := p[1:7];
  end for;
  // Standard errors from the curvature at the fit: the inverse of J'J + W in scaled coordinates.
  if iterations > 0 then
    J := zeros(nu * nv, np);
    for k in 1:np loop
      if free[k] then
        pTry := p;
        pTry[k] := p[k] + h * scale[k];
        for s in 1:ns loop
          srcTry[s, :] := pTry[7 + 7 * (s - 1) + 1:7 + 7 * (s - 1) + 7];
        end for;
        plus := renderImage(u, v, pTry[1:7], srcTry, pointMass, subpixels, psfSigma);
        pTry[k] := p[k] - h * scale[k];
        for s in 1:ns loop
          srcTry[s, :] := pTry[7 + 7 * (s - 1) + 1:7 + 7 * (s - 1) + 7];
        end for;
        minus := renderImage(u, v, pTry[1:7], srcTry, pointMass, subpixels, psfSigma);
        for i in 1:nu loop
          for j in 1:nv loop
            J[(i - 1) * nv + j, k] := -(plus[i, j] - minus[i, j]) / (2.0 * h * noiseSigma);
          end for;
        end for;
      end if;
    end for;
    A := transpose(J) * J;
    for k in 1:np loop
      A[k, k] := A[k, k] + weight[k];
      if not free[k] then
        A[k, :] := zeros(np);
        A[:, k] := zeros(np);
        A[k, k] := 1.0;
      end if;
    end for;
  end if;
  (C, ok) := LinearAlgebra.solve(A, identity(np));
  lens := p[1:7];
  for k in 1:7 loop
    lensSigma[k] := if free[k] and ok and iterations > 0 then sqrt(max(C[k, k], 0.0)) * scale[k] else 0.0;
  end for;
  for s in 1:ns loop
    for k in 1:7 loop
      sources[s, k] := p[7 + 7 * (s - 1) + k];
      sourcesSigma[s, k] := if free[7 + 7 * (s - 1) + k] and ok and iterations > 0 then sqrt(max(C[7 + 7 * (s - 1) + k, 7 + 7 * (s - 1) + k], 0.0)) * scale[7 + 7 * (s - 1) + k] else 0.0;
    end for;
    // Canonical source shape: axis ratio at most 1 (a larger one is the same ellipse with the axes swapped) and a position angle in [0, pi).
    if sources[s, 5] > 1.0 then
      sourcesSigma[s, 5] := sourcesSigma[s, 5] / sources[s, 5]^2;
      sources[s, 5] := 1.0 / sources[s, 5];
      sources[s, 6] := sources[s, 6] + 0.5 * GeneralRelativity.Constants.pi;
    end if;
    sources[s, 6] := sources[s, 6] - GeneralRelativity.Constants.pi * floor(sources[s, 6] / GeneralRelativity.Constants.pi);
  end for;
  // Canonical lens: a non-negative shear (a negative one is the same shear at 90 degrees) and position angles in [0, pi).
  if lens[6] < 0.0 then
    lens[6] := -lens[6];
    lens[7] := lens[7] + 0.5 * GeneralRelativity.Constants.pi;
  end if;
  lens[5] := lens[5] - GeneralRelativity.Constants.pi * floor(lens[5] / GeneralRelativity.Constants.pi);
  lens[7] := lens[7] - GeneralRelativity.Constants.pi * floor(lens[7] / GeneralRelativity.Constants.pi);
  annotation(Documentation(info="<html>
<p>Minimizes <code>&chi;&sup2; = &Sigma;<sub>pixels</sub> [(observed &minus; simulated)/&sigma;]&sup2; + &Sigma;<sub>k</sub> [(p<sub>k</sub> &minus; p<sub>k,0</sub>)/&sigma;<sub>k,prior</sub>]&sup2;</code>
over the free parameters. Each iteration builds the Jacobian of the pixel residuals by central differences (two simulations per
free parameter, in parameters scaled to order one), solves the damped normal equations
<code>(J<sup>T</sup>J + W + &lambda;D)&delta; = &minus;(J<sup>T</sup>r + W(p &minus; p<sub>0</sub>))</code> and accepts the step
when the cost falls, raising &lambda; tenfold otherwise (Marquardt's <code>D = diag(J<sup>T</sup>J)</code>).</p>
<p><b>Kalman view.</b> With a prior (<code>W = diag(1/&sigma;<sub>prior</sub>&sup2;)</code>) each accepted Gauss-Newton step is
the iterated extended Kalman filter's measurement update of a constant state (the lens) from one measurement vector (the pixels):
prior covariance <code>W<sup>&minus;1</sup></code>, measurement Jacobian <code>J</code>, unit measurement covariance in units of
&sigma;. The posterior covariance <code>(J<sup>T</sup>J + W)<sup>&minus;1</sup></code> gives the returned standard errors, so the
accuracy of the mass estimate is <code>2 &sigma;<sub>&theta;E</sub>/&theta;<sub>E</sub></code> relative (the mass goes as
&theta;<sub>E</sub>&sup2;). A prior is also how catalogue knowledge enters: a known lens position or mass is a tight prior on those
entries. Without priors the result is plain nonlinear least squares.</p>
<p>Parameters are kept physical after every step (&theta;<sub>E</sub> &ge; 0, 0.1 &le; q<sub>lens</sub> &le; 1, 0.1 &le; q<sub>source</sub> &le; 10 where q above 1 means the axes swapped, I<sub>e</sub> &ge; 0,
0.3 &le; n &le; 8). A point-mass lens fixes q and the position angle, which it does not use.</p>
</html>"));
end fitLens;
