within GeneralRelativity.Connections;
function curvatureAt "Connection coefficients, curvature, Ricci tensor and torsion of a connection function at one event"
  input GeneralRelativity.Interfaces.ConnectionFunction connection "Connection coefficients";
  input Real x[4] "Coordinates";
  input Real p[:] "Connection parameters";
  input Real step[4] "Finite-difference step per coordinate";
  output ConnectionCurvature k;
protected
  Real dGamma[4, 4, 4, 4];
algorithm
  k.x := x;
  k.Gamma := connection(x, p);
  dGamma := derivative(connection, x, p, step);
  k.riemann := riemann(k.Gamma, dGamma);
  k.ricci := ricci(k.riemann);
  k.torsion := torsion(k.Gamma);
end curvatureAt;
