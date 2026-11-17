within GeneralRelativity.Visualization;
impure function writeLensFit "Append the result of a lens inversion for the viewer's lensfit style: parameters before and after the fit, their standard errors, the truth when known, and the cost and mass history"
  input String fileName;
  input String lensModel = "galaxy" "\"galaxy\" (isothermal ellipse plus shear) or \"point\"";
  input String parameterNames[:] "One label per parameter";
  input Real initialValues[size(parameterNames, 1)] "Parameters before the fit (step 2 estimate), in display units";
  input Real fitted[size(parameterNames, 1)] "Parameters after the fit";
  input Real sigma[size(parameterNames, 1)] "Standard errors of the fitted parameters";
  input Real truth[:] = fill(0.0, 0) "True parameters when the image is synthetic; empty when unknown";
  input Real chi2History[:] "Cost before each iteration and after the last";
  input Real massHistorySolar[:] = fill(0.0, 0) "Projected mass inside the Einstein radius, in solar masses, per iteration";
  input Real massSolar = 0.0 "Projected mass inside the Einstein radius, fitted";
  input Real massSigmaSolar = 0.0 "Its standard error";
  input Real massInitialSolar = 0.0 "The step-2 estimate from the ring radius";
  input Real trueMassSolar = 0.0 "0 when unknown";
  input Real thetaEArcsec = 0.0 "Fitted Einstein radius";
  input Real distancesMpc[3] = zeros(3) "D_l, D_s, D_ls";
  input Integer iterations = 0;
  input Boolean converged = true;
  input Real pixels = 0.0 "Number of pixels, the expected chi2 of a correct model with the right noise";
  input Real explained = 0.0 "Fraction of the image explained by the fitted simulation (imageResidual)";
  input Real correlation = 0.0 "Correlation of observed and fitted images";
protected
  String names;
algorithm
  names := "[";
  for i in 1:size(parameterNames, 1) loop
    names := names + (if i > 1 then "," else "") + "\"" + parameterNames[i] + "\"";
  end for;
  names := names + "]";
  Modelica.Utilities.Streams.print("{\"type\":\"lensFit\",\"model\":\"" + lensModel + "\",\"parameterNames\":" + names
    + ",\"initial\":" + jsonVector(initialValues) + ",\"fitted\":" + jsonVector(fitted) + ",\"sigma\":" + jsonVector(sigma)
    + ",\"truth\":" + jsonVector(truth) + ",\"chi2History\":" + jsonVector(chi2History) + ",\"massHistorySolar\":" + jsonVector(massHistorySolar)
    + ",\"massSolar\":" + jsonNumber(massSolar) + ",\"massSigmaSolar\":" + jsonNumber(massSigmaSolar)
    + ",\"massInitialSolar\":" + jsonNumber(massInitialSolar) + ",\"trueMassSolar\":" + jsonNumber(trueMassSolar)
    + ",\"thetaEArcsec\":" + jsonNumber(thetaEArcsec) + ",\"distancesMpc\":" + jsonVector(distancesMpc)
    + ",\"iterations\":" + String(iterations) + ",\"converged\":" + (if converged then "true" else "false")
    + ",\"pixels\":" + jsonNumber(pixels) + ",\"explained\":" + jsonNumber(explained) + ",\"correlation\":" + jsonNumber(correlation) + "}", fileName);
end writeLensFit;
