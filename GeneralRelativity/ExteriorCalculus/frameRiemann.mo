within GeneralRelativity.ExteriorCalculus;
function frameRiemann "Riemann tensor with its first pair projected on a tetrad: R^a_b mu nu, the curvature 2-forms"
  input Real riemann[4, 4, 4, 4] "Covariant Riemann tensor R_lambda sigma mu nu";
  input Real e[4, 4] "Orthonormal tetrad e[a, mu] = e_a^mu";
  output Real Omega[4, 4, 4, 4] "Omega[a, b, mu, nu] = theta^a_rho R^rho_sigma mu nu e_b^sigma = eta^aa e_a^lambda e_b^sigma R_lambda sigma mu nu";
protected
  constant Real eta[4] = {-1.0, 1.0, 1.0, 1.0};
algorithm
  for a in 1:4 loop
    for b in 1:4 loop
      for mu in 1:4 loop
        for nu in 1:4 loop
          Omega[a, b, mu, nu] := 0.0;
          for l in 1:4 loop
            for s in 1:4 loop
              Omega[a, b, mu, nu] := Omega[a, b, mu, nu] + eta[a] * e[a, l] * e[b, s] * riemann[l, s, mu, nu];
            end for;
          end for;
        end for;
      end for;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>For each frame pair <code>(a, b)</code>, <code>Omega[a, b, :, :]</code> is the
coordinate-component 2-form <code>&Omega;<sup>a</sup><sub>b</sub> = &frac12;R<sup>a</sup><sub>b&mu;&nu;</sub>dx<sup>&mu;</sup>&and;dx<sup>&nu;</sup></code>
that Cartan's second structure equation produces
(<a href=\"modelica://GeneralRelativity.ExteriorCalculus.curvatureForm\">curvatureForm</a>).</p>
</html>"));
end frameRiemann;
