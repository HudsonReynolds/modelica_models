within GeneralRelativity.Lensing.Inversion;
function deflectionAt "Scaled deflection in the sky frame of a lens parameter vector: point mass or isothermal ellipse, plus external shear"
  input Real theta[2](each unit = "rad") "Image-plane position";
  input Real lens[7] "thetaE, x0, y0, q, position angle, shear, shear angle (Inversion package documentation)";
  input Boolean pointMass = false "true: point mass (q and the position angle are ignored)";
  output Real alpha[2](each unit = "rad");
protected
  Real d[2] = theta - lens[2:3];
  Real c;
  Real s;
  Real aw[2];
algorithm
  if pointMass then
    alpha := GeneralRelativity.Lensing.deflectionPointMass(d, abs(lens[1]));
  else
    c := cos(lens[5]);
    s := sin(lens[5]);
    aw := GeneralRelativity.Lensing.deflectionSIE({c * d[1] + s * d[2], -s * d[1] + c * d[2]}, abs(lens[1]), min(max(lens[4], 1.0e-3), 1.0));
    alpha := {c * aw[1] - s * aw[2], s * aw[1] + c * aw[2]};
  end if;
  alpha := alpha + GeneralRelativity.Lensing.deflectionShear(d, lens[6], lens[7]);
end deflectionAt;
