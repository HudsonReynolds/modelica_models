within GeneralRelativity.Tensors;
function ricciSquared "Quadratic Ricci invariant R_ab R^ab"
  input Real gInv[:, size(gInv, 1)] "Contravariant metric";
  input Real Ric[size(gInv, 1), size(gInv, 1)] "Ricci tensor R_ab";
  output Real Ric2 "R_ab R^ab";
protected
  Integer n = size(gInv, 1);
  Real mixed[size(gInv, 1), size(gInv, 1)] "R^a_b";
algorithm
  mixed := gInv * Ric;
  Ric2 := sum(mixed[a, b] * mixed[b, a] for a in 1:n, b in 1:n);
  annotation(Documentation(info="<html>
<p>Computed as the trace of the square of the mixed tensor <code>R<sup>a</sup><sub>b</sub></code>.
It is not sign-definite in Lorentzian signature.</p>
</html>"));
end ricciSquared;
