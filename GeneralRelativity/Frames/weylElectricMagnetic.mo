within GeneralRelativity.Frames;
function weylElectricMagnetic "Electric and magnetic parts E_ij = C_0i0j, B_ij = (1/2) eps_ikl C_kl0j of the Weyl tensor in the Eulerian orthonormal frame"
  input Real g[4, 4] "Covariant metric";
  input Real C[4, 4, 4, 4] "Weyl tensor C_abcd (coordinate components)";
  output Real E[3, 3] "Electric part E_ij (symmetric, traceless)";
  output Real B[3, 3] "Magnetic part B_ij (symmetric, traceless)";
protected
  Real F[4, 4, 4, 4];
algorithm
  F := frameComponents4(C, orthonormalTetrad(g));
  for i in 1:3 loop
    for j in 1:3 loop
      E[i, j] := F[1, i + 1, 1, j + 1];
      // Levi-Civita symbol eps_ikl = (i - k)(k - l)(l - i)/2 for i, k, l in 1..3.
      B[i, j] := 0.5 * sum((i - k) * (k - m) * (m - i) / 2 * F[k + 1, m + 1, 1, j + 1] for k in 1:3, m in 1:3);
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>Frame indices refer to <a href=\"modelica://GeneralRelativity.Frames.orthonormalTetrad\">orthonormalTetrad</a>
(index 0 is the Eulerian observer). <code>E<sub>ij</sub></code> is the tidal tensor felt by that observer,
<code>D&sup2;&xi;<sup>i</sup>/d&tau;&sup2; = &minus;E<sup>i</sup><sub>j</sub>&xi;<sup>j</sup></code> in vacuum;
<code>B<sub>ij</sub></code> is the frame-dragging (gravitomagnetic) part. In four dimensions they determine
the Weyl tensor completely, and the complex matrix <code>Q = E + iB</code> carries its Petrov type
(<a href=\"modelica://GeneralRelativity.Frames.classifyPetrovMatrix\">classifyPetrovMatrix</a>).</p>
</html>"));
end weylElectricMagnetic;
