within GeneralRelativity.Frames;
function frameComponents2 "Frame components T_ab = e_a^mu e_b^nu T_mu nu of a covariant rank-2 tensor"
  input Real T[4, 4] "Covariant tensor T_mu nu";
  input Real e[4, 4] "Tetrad, e[a, :] = e_a^mu";
  output Real Tf[4, 4] "Frame components T_ab";
algorithm
  Tf := e * T * transpose(e);
  annotation(Documentation(info="<html>
<p>With an orthonormal tetrad, <code>frameComponents2(g, e)</code> is <code>&eta;</code> and the
frame components of a stress-energy tensor are the energy density <code>T<sub>00</sub></code>,
momentum density <code>&minus;T<sub>0i</sub></code> and stresses <code>T<sub>ij</sub></code> measured by observer
<code>e<sub>0</sub></code>.</p>
</html>"));
end frameComponents2;
