within GeneralRelativity.Perturbation;
function traceReverse "Trace-reversed perturbation hbar_ab = h_ab - 1/2 eta_ab h, an involution"
  input Real h[4, 4] "Metric perturbation h_ab";
  output Real hBar[4, 4];
protected
  constant Real eta[4, 4] = [-1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
  Real traceH = sum(eta[a, a] * h[a, a] for a in 1:4);
algorithm
  hBar := h - 0.5 * traceH * eta;
end traceReverse;
