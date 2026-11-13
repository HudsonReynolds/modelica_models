within GeneralRelativity.Tensors;
function fullContraction "Full contraction A_abcd B^abcd of two covariant rank-4 tensors, in any dimension"
  input Real gInv[:, size(gInv, 1)] "Contravariant metric";
  input Real A[size(gInv, 1), size(gInv, 1), size(gInv, 1), size(gInv, 1)] "A_abcd";
  input Real B[size(gInv, 1), size(gInv, 1), size(gInv, 1), size(gInv, 1)] "B_abcd";
  output Real AB "A_abcd B^abcd";
protected
  Integer n = size(gInv, 1);
  Real B1[size(gInv, 1), size(gInv, 1), size(gInv, 1), size(gInv, 1)];
  Real B2[size(gInv, 1), size(gInv, 1), size(gInv, 1), size(gInv, 1)];
algorithm
  // Raise one index at a time: four passes of n^5 products instead of n^8.
  for a in 1:n loop
    for b in 1:n loop
      for c in 1:n loop
        for d in 1:n loop
          B1[a, b, c, d] := sum(gInv[a, e] * B[e, b, c, d] for e in 1:n);
        end for;
      end for;
    end for;
  end for;
  for a in 1:n loop
    for b in 1:n loop
      for c in 1:n loop
        for d in 1:n loop
          B2[a, b, c, d] := sum(gInv[b, e] * B1[a, e, c, d] for e in 1:n);
        end for;
      end for;
    end for;
  end for;
  for a in 1:n loop
    for b in 1:n loop
      for c in 1:n loop
        for d in 1:n loop
          B1[a, b, c, d] := sum(gInv[c, e] * B2[a, b, e, d] for e in 1:n);
        end for;
      end for;
    end for;
  end for;
  for a in 1:n loop
    for b in 1:n loop
      for c in 1:n loop
        for d in 1:n loop
          B2[a, b, c, d] := sum(gInv[d, e] * B1[a, b, c, e] for e in 1:n);
        end for;
      end for;
    end for;
  end for;
  AB := sum(A[a, b, c, d] * B2[a, b, c, d] for a in 1:n, b in 1:n, c in 1:n, d in 1:n);
  annotation(Documentation(info="<html>
<p><code>fullContraction(gInv, R, R)</code> is the
<a href=\"modelica://GeneralRelativity.Tensors.kretschmann\">Kretschmann scalar</a> and
<code>fullContraction(gInv, C, C)</code> with the <a href=\"modelica://GeneralRelativity.Tensors.weyl\">Weyl tensor</a>
is the Weyl invariant <code>C<sub>abcd</sub>C<sup>abcd</sup></code>.</p>
</html>"));
end fullContraction;
