within GeneralRelativity.Visualization;
impure function writeLens "Declare a gravitational-lens scene for the viewer's lensing style: lens model, mass, redshifts, cosmology and the Einstein radius computed in Modelica"
  input String fileName;
  input String lensModel = "galaxy" "\"galaxy\" (singular isothermal ellipse plus shear) or \"point\" (point mass)";
  input Real massSolar = 5.0e11 "Lens mass in solar masses: total for a point, projected within the aperture for a galaxy";
  input Real apertureKpc = 10.0 "Aperture radius that defines the galaxy mass";
  input Real zLens = 0.5;
  input Real zSource = 2.0;
  input Real H0 = 70.0 "km/s/Mpc";
  input Real OmegaM = 0.3 "Flat universe: OmegaLambda = 1 - OmegaM";
  input Real axisRatio = 0.75 "Galaxy axis ratio q";
  input Real positionAngleDeg = 30.0 "Galaxy major-axis position angle";
  input Real shear = 0.05 "External shear";
  input Real shearAngleDeg = 110.0 "External shear position angle";
  input Real thetaEArcsec = 0.0 "Einstein radius computed in Modelica, shown by the viewer for comparison";
  input Real distancesMpc[3] = zeros(3) "D_l, D_s, D_ls computed in Modelica";
algorithm
  Modelica.Utilities.Streams.print("{\"type\":\"lens\",\"model\":\"" + lensModel + "\",\"massSolar\":" + jsonNumber(massSolar)
    + ",\"apertureKpc\":" + jsonNumber(apertureKpc) + ",\"zLens\":" + jsonNumber(zLens) + ",\"zSource\":" + jsonNumber(zSource)
    + ",\"H0\":" + jsonNumber(H0) + ",\"OmegaM\":" + jsonNumber(OmegaM) + ",\"axisRatio\":" + jsonNumber(axisRatio)
    + ",\"positionAngleDeg\":" + jsonNumber(positionAngleDeg) + ",\"shear\":" + jsonNumber(shear)
    + ",\"shearAngleDeg\":" + jsonNumber(shearAngleDeg) + ",\"thetaEArcsec\":" + jsonNumber(thetaEArcsec)
    + ",\"distancesMpc\":" + jsonVector(distancesMpc) + "}", fileName);
end writeLens;
