within GeneralRelativity.Connections;
record ConnectionCurvature "Coefficients, curvature and torsion of an affine connection at one event"
  Real x[4] "Coordinates of the event";
  Real Gamma[4, 4, 4] "Gamma^i_jk";
  Real riemann[4, 4, 4, 4] "R^r_smn";
  Real ricci[4, 4] "R_sn = R^r_srn";
  Real torsion[4, 4, 4] "T^i_jk";
end ConnectionCurvature;
