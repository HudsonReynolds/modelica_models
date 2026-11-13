within GeneralRelativity.Tensors;
function kretschmann "Kretschmann invariant K = R_abcd R^abcd"
  input Real gInv[:, size(gInv, 1)] "Contravariant metric";
  input Real RLower[size(gInv, 1), size(gInv, 1), size(gInv, 1), size(gInv, 1)] "R_abcd";
  output Real K "Kretschmann scalar";
protected
  Integer n = size(gInv, 1);
  Real R1[size(gInv, 1), size(gInv, 1), size(gInv, 1), size(gInv, 1)];
  Real R2[size(gInv, 1), size(gInv, 1), size(gInv, 1), size(gInv, 1)];
algorithm
  // Raise one index at a time: four passes of n^5 products instead of n^8.
  for a in 1:n loop
    for b in 1:n loop
      for c in 1:n loop
        for d in 1:n loop
          R1[a, b, c, d] := sum(gInv[a, e] * RLower[e, b, c, d] for e in 1:n);
        end for;
      end for;
    end for;
  end for;
  for a in 1:n loop
    for b in 1:n loop
      for c in 1:n loop
        for d in 1:n loop
          R2[a, b, c, d] := sum(gInv[b, e] * R1[a, e, c, d] for e in 1:n);
        end for;
      end for;
    end for;
  end for;
  for a in 1:n loop
    for b in 1:n loop
      for c in 1:n loop
        for d in 1:n loop
          R1[a, b, c, d] := sum(gInv[c, e] * R2[a, b, e, d] for e in 1:n);
        end for;
      end for;
    end for;
  end for;
  for a in 1:n loop
    for b in 1:n loop
      for c in 1:n loop
        for d in 1:n loop
          R2[a, b, c, d] := sum(gInv[d, e] * R1[a, b, c, e] for e in 1:n);
        end for;
      end for;
    end for;
  end for;
  K := sum(RLower[a, b, c, d] * R2[a, b, c, d] for a in 1:n, b in 1:n, c in 1:n, d in 1:n);
  annotation(Documentation(info="<html>
<p>A coordinate-independent measure of curvature that does not vanish in vacuum.
For Schwarzschild <code>K = 12 r<sub>s</sub>&sup2;/r<sup>6</sup> = 48 M&sup2;/r<sup>6</sup></code>,
finite at the horizon and divergent at <code>r = 0</code>; the visualization
package colours embedding surfaces with it.</p>
</html>"));
end kretschmann;
