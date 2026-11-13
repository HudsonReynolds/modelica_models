within GeneralRelativity.Tensors;
function curvatureInvariants "Quadratic and cubic curvature invariants of a Curvature record"
  input Curvature k "Curvature at one event, e.g. from Metrics.curvatureAt";
  output Invariants inv "Scalar invariants";
protected
  Real C[4, 4, 4, 4] "Weyl tensor";
algorithm
  C := weyl(k.g, k.riemann, k.ricci, k.ricciScalar);
  inv.ricciScalar := k.ricciScalar;
  inv.ricciSquared := ricciSquared(k.gInv, k.ricci);
  inv.kretschmann := k.kretschmann;
  inv.weylSquared := fullContraction(k.gInv, C, C);
  inv.gaussBonnet := gaussBonnet(inv.kretschmann, inv.ricciSquared, inv.ricciScalar);
  inv.cubic := cubicInvariants(k.gInv, k.riemann, k.ricci, k.ricciScalar);
  inv.weylCubed := weylCubed(k.gInv, C);
end curvatureInvariants;
