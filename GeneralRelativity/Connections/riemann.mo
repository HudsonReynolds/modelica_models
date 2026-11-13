within GeneralRelativity.Connections;
function riemann "Curvature R[r, s, m, n] = R^r_smn of any affine connection, with or without torsion"
  input Real Gamma[4, 4, 4] "Gamma[i, j, k] = Gamma^i_jk, derivative index j";
  input Real dGamma[4, 4, 4, 4] "dGamma[m, i, j, k] = d_m Gamma^i_jk";
  output Real R[4, 4, 4, 4] "R^r_smn, antisymmetric in m, n";
algorithm
  for r in 1:4 loop
    for s in 1:4 loop
      for m in 1:4 loop
        for n in 1:4 loop
          R[r, s, m, n] := dGamma[m, r, n, s] - dGamma[n, r, m, s]
            + sum(Gamma[r, m, l] * Gamma[l, n, s] - Gamma[r, n, l] * Gamma[l, m, s] for l in 1:4);
        end for;
      end for;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p><code>R<sup>&rho;</sup><sub>&sigma;&mu;&nu;</sub> = &part;<sub>&mu;</sub>&Gamma;<sup>&rho;</sup><sub>&nu;&sigma;</sub> &minus; &part;<sub>&nu;</sub>&Gamma;<sup>&rho;</sup><sub>&mu;&sigma;</sub> + &Gamma;<sup>&rho;</sup><sub>&mu;&lambda;</sub>&Gamma;<sup>&lambda;</sup><sub>&nu;&sigma;</sub> &minus; &Gamma;<sup>&rho;</sup><sub>&nu;&lambda;</sub>&Gamma;<sup>&lambda;</sup><sub>&mu;&sigma;</sub></code>.
For the Levi-Civita connection it equals
<a href=\"modelica://GeneralRelativity.Tensors.riemannMixed\">Tensors.riemannMixed</a>
of <a href=\"modelica://GeneralRelativity.Tensors.riemannLower\">Tensors.riemannLower</a>.
It needs first derivatives of the connection, so a metric connection built by
finite differences needs nested differences of the metric; prefer
<a href=\"modelica://GeneralRelativity.Metrics.curvatureAt\">Metrics.curvatureAt</a> there.</p>
</html>"));
end riemann;
