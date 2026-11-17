within GeneralRelativity.Cosmology;
function reionizationRedshift "Reionization redshift that gives the optical depth par.tauReio (bisection)"
  input Parameters par;
  output Real zReio;
protected
  Real lo = 2.0;
  Real hi = 30.0;
  Real tau;
algorithm
  for iter in 1:60 loop
    zReio := 0.5 * (lo + hi);
    tau := reionizationOpticalDepth(zReio, par);
    if tau > par.tauReio then
      hi := zReio;
    else
      lo := zReio;
    end if;
  end for;
  zReio := 0.5 * (lo + hi);
end reionizationRedshift;
