within GeneralRelativity.Cosmology;
impure function writeExpansionRecords "Write the cosmology and expansion-history records of the expanding-universe scene; returns the key epochs"
  input Parameters par;
  input String fileName;
  input Real zStar "Redshift of last scattering (from thermalHistory)";
  input Integer n = 400 "Rows, logarithmic in time";
  input Real tMin(unit = "Gyr") = 1.0e-7;
  input Real tMax(unit = "Gyr") = 60.0;
  output Real epochs[12] "age, zEq, tEq, tStar, zLambda, tLambda, zAcc, tAcc, chiStar (Gly), hubbleRadiusToday, particleHorizonToday, eventHorizonToday (Gly)";
protected
  Real H0, Og, On, Ob, Oc, OL;
  Real t[n], a[n], z[n], H[n], T[n], eta[n], chi[n], dH[n], dP[n], dE[n], q[n];
  Real eta0, age, zEq, zLambda, zAcc, lo, hi, mid;
algorithm
  (H0, Og, On, Ob, Oc, OL) := densities(par);
  age := cosmicTime(1.0, par);
  eta0 := conformalTime(1.0, par) * MpcGyr;
  for i in 1:n loop
    t[i] := tMin * (tMax / tMin)^((i - 1) / (n - 1));
    a[i] := scaleFactorAtTime(t[i], par);
    z[i] := 1.0 / a[i] - 1.0;
    H[i] := hubble(a[i], par);
    T[i] := par.TCMB / a[i];
    eta[i] := conformalTime(a[i], par) * MpcGyr;
    chi[i] := eta0 - eta[i];
    dH[i] := GeneralRelativity.Constants.c / 1000.0 / H[i] * MpcGyr;
    dP[i] := a[i] * eta[i];
    dE[i] := a[i] * eventHorizon(a[i], par) * MpcGyr;
    q[i] := decelerationParameter(a[i], par);
  end for;
  zEq := (Ob + Oc) / (Og + On) - 1.0;
  zLambda := (OL / (Ob + Oc))^(1.0 / 3.0) - 1.0;
  // The expansion starts to accelerate where q = 0.
  lo := 0.1;
  hi := 1.0;
  for i in 1:60 loop
    mid := 0.5 * (lo + hi);
    if decelerationParameter(mid, par) > 0.0 then
      lo := mid;
    else
      hi := mid;
    end if;
  end for;
  zAcc := 1.0 / (0.5 * (lo + hi)) - 1.0;
  epochs := {age, zEq, cosmicTime(1.0 / (1.0 + zEq), par), cosmicTime(1.0 / (1.0 + zStar), par), zLambda, cosmicTime(1.0 / (1.0 + zLambda), par),
    zAcc, cosmicTime(1.0 / (1.0 + zAcc), par), eta0 - conformalTime(1.0 / (1.0 + zStar), par) * MpcGyr,
    GeneralRelativity.Constants.c / 1000.0 / hubble(1.0, par) * MpcGyr, eta0, eventHorizon(1.0, par) * MpcGyr};
  Modelica.Utilities.Streams.print("{\"type\":\"cosmology\",\"h\":" + String(par.h) + ",\"omegaB\":" + String(Ob) + ",\"omegaC\":" + String(Oc)
    + ",\"omegaM\":" + String(Ob + Oc) + ",\"omegaR\":" + String(Og + On) + ",\"omegaL\":" + String(OL) + ",\"omegaK\":" + String(par.omegaK)
    + ",\"Tcmb\":" + String(par.TCMB) + ",\"age\":" + String(age) + ",\"zEq\":" + String(zEq) + ",\"tEq\":" + String(epochs[3])
    + ",\"zStar\":" + String(zStar) + ",\"tStar\":" + String(epochs[4]) + ",\"zLambda\":" + String(zLambda) + ",\"tLambda\":" + String(epochs[6])
    + ",\"zAcc\":" + String(zAcc) + ",\"tAcc\":" + String(epochs[8]) + ",\"chiStar\":" + String(epochs[9])
    + ",\"hubbleRadiusToday\":" + String(epochs[10]) + ",\"particleHorizonToday\":" + String(epochs[11])
    + ",\"eventHorizonToday\":" + String(epochs[12]) + "}", fileName);
  Modelica.Utilities.Streams.print("{\"type\":\"expansion\",\"t\":" + GeneralRelativity.Visualization.jsonVector(t, 7)
    + ",\"a\":" + GeneralRelativity.Visualization.jsonVector(a, 7) + ",\"z\":" + GeneralRelativity.Visualization.jsonVector(z, 7)
    + ",\"H\":" + GeneralRelativity.Visualization.jsonVector(H, 7) + ",\"T\":" + GeneralRelativity.Visualization.jsonVector(T, 7)
    + ",\"eta\":" + GeneralRelativity.Visualization.jsonVector(eta, 7) + ",\"chi\":" + GeneralRelativity.Visualization.jsonVector(chi, 7)
    + ",\"dHubble\":" + GeneralRelativity.Visualization.jsonVector(dH, 7) + ",\"dParticle\":" + GeneralRelativity.Visualization.jsonVector(dP, 7)
    + ",\"dEvent\":" + GeneralRelativity.Visualization.jsonVector(dE, 7) + ",\"q\":" + GeneralRelativity.Visualization.jsonVector(q, 7) + "}", fileName);
  annotation(Documentation(info="<html>
<p>Rows are logarithmic in cosmic time <code>t</code> (Gyr). Distances are in Gly (light travels 1 Gly per Gyr):
<code>eta</code> conformal time (comoving distance light has covered since the big bang), <code>chi</code> the comoving
distance to a source whose light emitted at <code>t</code> arrives today (negative in the future),
<code>dHubble</code> the proper Hubble radius <code>c/H</code>, <code>dParticle</code> the proper particle horizon
<code>a&eta;</code> and <code>dEvent</code> the proper event horizon.</p>
</html>"));
end writeExpansionRecords;
