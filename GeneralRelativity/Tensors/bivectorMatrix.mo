within GeneralRelativity.Tensors;
function bivectorMatrix "Rank-4 tensor C_abcd as the matrix M[(a-1) n + b, (c-1) n + d] = C_ab^cd acting on index pairs"
  input Real gInv[:, size(gInv, 1)] "Contravariant metric";
  input Real C[size(gInv, 1), size(gInv, 1), size(gInv, 1), size(gInv, 1)] "C_abcd";
  output Real M[size(gInv, 1) * size(gInv, 1), size(gInv, 1) * size(gInv, 1)] "C_ab^cd";
protected
  Integer n = size(gInv, 1);
  Real C1[size(gInv, 1), size(gInv, 1), size(gInv, 1), size(gInv, 1)] "C_ab^c_d";
algorithm
  for a in 1:n loop
    for b in 1:n loop
      for c in 1:n loop
        for d in 1:n loop
          C1[a, b, c, d] := sum(gInv[c, e] * C[a, b, e, d] for e in 1:n);
        end for;
      end for;
    end for;
  end for;
  for a in 1:n loop
    for b in 1:n loop
      for c in 1:n loop
        for d in 1:n loop
          M[(a - 1) * n + b, (c - 1) * n + d] := sum(gInv[d, e] * C1[a, b, c, e] for e in 1:n);
        end for;
      end for;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>Traces of powers of <code>M</code> give the invariants built by chaining index pairs:
<code>tr M&sup2; = C<sub>abcd</sub>C<sup>cdab</sup></code> and
<code>tr M&sup3; = C<sub>ab</sub><sup>cd</sup>C<sub>cd</sub><sup>ef</sup>C<sub>ef</sub><sup>ab</sup></code>.
For a tensor antisymmetric in each pair, <code>M</code> annihilates symmetric pairs and acts on a
bivector <code>F<sub>cd</sub></code> as <code>F<sub>ab</sub> &rarr; C<sub>ab</sub><sup>cd</sup>F<sub>cd</sub></code>.</p>
</html>"));
end bivectorMatrix;
