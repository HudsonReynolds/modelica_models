within GeneralRelativity.Tensors;
function traceFreeRicci "Trace-free Ricci tensor S_ab = R_ab - R g_ab / n"
  input Real g[:, size(g, 1)] "Covariant metric";
  input Real Ric[size(g, 1), size(g, 1)] "Ricci tensor R_ab";
  input Real R "Ricci scalar";
  output Real S[size(g, 1), size(g, 1)] "S_ab";
algorithm
  S := Ric - R / size(g, 1) * g;
  annotation(Documentation(info="<html>
<p>Zero exactly for Einstein spaces <code>R<sub>ab</sub> = &Lambda;g<sub>ab</sub></code>
(vacuum with or without a cosmological constant). In four dimensions it is proportional to the
trace-free part of the stress-energy, <code>S<sub>ab</sub> = 8&pi;G/c<sup>4</sup> (T<sub>ab</sub> &minus; T g<sub>ab</sub>/4)</code>.</p>
</html>"));
end traceFreeRicci;
