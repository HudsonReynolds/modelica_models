within GeneralRelativity.Frames;
function frameComponents4 "Frame components T_abcd = e_a^mu e_b^nu e_c^rho e_d^sigma T_mu nu rho sigma of a covariant rank-4 tensor"
  input Real T[4, 4, 4, 4] "Covariant tensor, e.g. R_abcd";
  input Real e[4, 4] "Tetrad, e[a, :] = e_a^mu";
  output Real Tf[4, 4, 4, 4] "Frame components";
protected
  Real T1[4, 4, 4, 4];
algorithm
  // Project one index at a time: four passes of 4^5 products.
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        for d in 1:4 loop
          T1[a, b, c, d] := sum(e[a, m] * T[m, b, c, d] for m in 1:4);
        end for;
      end for;
    end for;
  end for;
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        for d in 1:4 loop
          Tf[a, b, c, d] := sum(e[b, m] * T1[a, m, c, d] for m in 1:4);
        end for;
      end for;
    end for;
  end for;
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        for d in 1:4 loop
          T1[a, b, c, d] := sum(e[c, m] * Tf[a, b, m, d] for m in 1:4);
        end for;
      end for;
    end for;
  end for;
  for a in 1:4 loop
    for b in 1:4 loop
      for c in 1:4 loop
        for d in 1:4 loop
          Tf[a, b, c, d] := sum(e[d, m] * T1[a, b, c, m] for m in 1:4);
        end for;
      end for;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>The rows of <code>e</code> need not be orthonormal: any four vectors may be used, e.g. the real
null basis <code>{l, n, Re m, Im m}</code>. For the Schwarzschild Riemann tensor in the
<a href=\"modelica://GeneralRelativity.Frames.orthonormalTetrad\">static tetrad</a>,
<code>R<sub>0101</sub> = &minus;r<sub>s</sub>/r&sup3;</code> (radial stretching).</p>
</html>"));
end frameComponents4;
