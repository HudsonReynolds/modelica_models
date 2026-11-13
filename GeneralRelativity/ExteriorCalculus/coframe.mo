within GeneralRelativity.ExteriorCalculus;
function coframe "Dual coframe theta^a_mu = eta_ab g_mu nu e_b^nu of an orthonormal tetrad"
  input Real g[4, 4] "Covariant metric";
  input Real e[4, 4] "Orthonormal tetrad e[a, mu] = e_a^mu";
  output Real theta[4, 4] "theta[a, mu] = theta^a_mu; theta^a(e_b) = delta^a_b";
protected
  constant Real eta[4] = {-1.0, 1.0, 1.0, 1.0};
algorithm
  for a in 1:4 loop
    theta[a, :] := eta[a] * (g * e[a, :]);
  end for;
  annotation(Documentation(info="<html>
<p>Lowers the frame vectors with the metric and the frame index with
<code>&eta;</code>, which avoids a matrix inverse; it equals the inverse of
<code>e</code> only for an orthonormal tetrad (check with
<a href=\"modelica://GeneralRelativity.ExteriorCalculus.orthonormalityResidual\">orthonormalityResidual</a>).
Then <code>g<sub>&mu;&nu;</sub> = &eta;<sub>ab</sub>&theta;<sup>a</sup><sub>&mu;</sub>&theta;<sup>b</sup><sub>&nu;</sub></code>.</p>
</html>"));
end coframe;
