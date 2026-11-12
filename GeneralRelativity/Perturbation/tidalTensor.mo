within GeneralRelativity.Perturbation;
function tidalTensor "Electric part of the Riemann tensor for a static observer, E_ij = c^2 R_0i0j; the Newtonian tidal tensor Phi_,ij"
  input Real R[4, 4, 4, 4] "R_abcd (full or linearized)";
  input Real c(unit = "m/s") = 1.0 "Speed of light";
  output Real E[3, 3](each unit = "1/s2") "Relative acceleration per unit separation, d^2 xi_i / dt^2 = -E_ij xi_j";
algorithm
  for i in 1:3 loop
    for j in 1:3 loop
      E[i, j] := c^2 * R[1, i + 1, 1, j + 1];
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>In a coordinate chart where <code>x[1] = ct</code> and the observer is at rest
with <code>u = (c, 0, 0, 0)</code>, the geodesic deviation equation reads
<code>d&sup2;&xi;<sup>i</sup>/d&tau;&sup2; = -c&sup2;R<sup>i</sup><sub>0j0</sub>&xi;<sup>j</sup></code>.
For a Newtonian point mass the eigenvalues are
<code>-2GM/r&sup3;</code> radially and <code>+GM/r&sup3;</code> transversely, traceless in
vacuum.</p>
</html>"));
end tidalTensor;
