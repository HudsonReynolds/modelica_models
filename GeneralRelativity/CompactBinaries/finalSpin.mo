within GeneralRelativity.CompactBinaries;
function finalSpin "Dimensionless spin of the black hole left by the merger of two non-spinning black holes"
  input Real eta "Symmetric mass ratio";
  output Real chi "a/M of the remnant";
algorithm
  chi := 2.0 * sqrt(3.0) * eta - 3.5171 * eta^2 + 2.5763 * eta^3;
  annotation(Documentation(info="<html>
<p>Fit to numerical-relativity simulations, L. Rezzolla et al., <i>Astrophys. J. Lett.</i> 674, L29 (2008),
for non-spinning progenitors: <code>&chi; = 2&radic;3 &eta; &minus; 3.5171&eta;&sup2; + 2.5763&eta;&sup3;</code>;
equal masses give 0.686. The leading term is the orbital angular momentum at the innermost stable orbit of
a test particle.</p>
</html>"));
end finalSpin;
