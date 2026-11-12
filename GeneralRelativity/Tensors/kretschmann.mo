within GeneralRelativity.Tensors;
function kretschmann "Kretschmann invariant K = R_abcd R^abcd"
  input Real gInv[4, 4] "Contravariant metric";
  input Real RLower[4, 4, 4, 4] "R_abcd";
  output Real K "Kretschmann scalar";
protected
  Real R1[4, 4, 4, 4];
  Real R2[4, 4, 4, 4];
algorithm
  // Raise one index at a time: four passes of 4^5 products instead of 4^8.
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        for d in 1:4 loop
          R1[a, b, c, d] := sum(gInv[a, e] * RLower[e, b, c, d] for e in 1:4);
        end for;
      end for;
    end for;
  end for;
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        for d in 1:4 loop
          R2[a, b, c, d] := sum(gInv[b, e] * R1[a, e, c, d] for e in 1:4);
        end for;
      end for;
    end for;
  end for;
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        for d in 1:4 loop
          R1[a, b, c, d] := sum(gInv[c, e] * R2[a, b, e, d] for e in 1:4);
        end for;
      end for;
    end for;
  end for;
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        for d in 1:4 loop
          R2[a, b, c, d] := sum(gInv[d, e] * R1[a, b, c, e] for e in 1:4);
        end for;
      end for;
    end for;
  end for;
  K := sum(RLower[a, b, c, d] * R2[a, b, c, d] for a in 1:4, b in 1:4, c in 1:4, d in 1:4);
  annotation(Documentation(info="<html>
<p>A coordinate-independent measure of curvature that does not vanish in vacuum.
For Schwarzschild <code>K = 12 r<sub>s</sub>&sup2;/r<sup>6</sup> = 48 M&sup2;/r<sup>6</sup></code>,
finite at the horizon and divergent at <code>r = 0</code>; the visualization
package colours embedding surfaces with it.</p>
</html>"));
end kretschmann;
