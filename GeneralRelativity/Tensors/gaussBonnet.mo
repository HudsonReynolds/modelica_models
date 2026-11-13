within GeneralRelativity.Tensors;
function gaussBonnet "Gauss-Bonnet (Euler) density K - 4 R_ab R^ab + R^2"
  input Real K "Kretschmann scalar R_abcd R^abcd";
  input Real Ric2 "Ricci invariant R_ab R^ab";
  input Real R "Ricci scalar";
  output Real E "Gauss-Bonnet density";
algorithm
  E := K - 4.0 * Ric2 + R^2;
  annotation(Documentation(info="<html>
<p>In four dimensions <code>&int;E &radic;(&minus;g) d<sup>4</sup>x</code> is topological
(<code>32&pi;&sup2;&chi;</code> for a closed Riemannian manifold) and its variation vanishes;
equivalently <code>E = C<sub>abcd</sub>C<sup>abcd</sup> &minus; 2R<sub>ab</sub>R<sup>ab</sup> + 2R&sup2;/3</code>.
In vacuum it equals the Kretschmann scalar.</p>
</html>"));
end gaussBonnet;
