within GeneralRelativity.IndexNotation;
function systemFromCurvature "System with the metric, Christoffel symbols, Riemann, Ricci, Einstein, Weyl, Schouten, trace-free Ricci and stress-energy tensors, the Ricci and Kretschmann scalars and epsilon of a Curvature record"
  input GeneralRelativity.Tensors.Curvature k "Curvature at one event, from Metrics.curvatureAt";
  input Real gravitationalConstant(unit = "m3/(kg.s2)") = 1.0 "Newton's constant for T; 1 in geometric units";
  input Real c(unit = "m/s") = 1.0 "Speed of light for T; 1 in geometric units";
  input Real Lambda(unit = "1/m2") = 0.0 "Cosmological constant for T";
  output TensorSystem sys;
protected
  Real GammaLower[4, 4, 4] "g_il Gamma^l_jk";
  Real C[4, 4, 4, 4];
algorithm
  sys := systemFromMetric(k.g);
  sys.gInv := k.gInv;
  for i in 1:4 loop
    for j in 1:4 loop
      for l in 1:4 loop
        GammaLower[i, j, l] := sum(k.g[i, m] * k.Gamma[m, j, l] for m in 1:4);
      end for;
    end for;
  end for;
  sys := Internal.register(sys, "Gamma", 3, Internal.flatten3(GammaLower));
  sys := Internal.register(sys, "R", 4, Internal.flatten4(k.riemann));
  sys := Internal.alias(sys, "Riemann", "R", 4);
  sys := Internal.register(sys, "R", 2, Internal.flatten2(k.ricci));
  sys := Internal.alias(sys, "Ric", "R", 2);
  sys := Internal.alias(sys, "Ricci", "R", 2);
  sys := Internal.register(sys, "R", 0, {k.ricciScalar});
  sys := Internal.register(sys, "G", 2, Internal.flatten2(k.einstein));
  sys := Internal.alias(sys, "Einstein", "G", 2);
  C := GeneralRelativity.Tensors.weyl(k.g, k.riemann, k.ricci, k.ricciScalar);
  sys := Internal.register(sys, "C", 4, Internal.flatten4(C));
  sys := Internal.alias(sys, "Weyl", "C", 4);
  sys := Internal.register(sys, "P", 2, Internal.flatten2(GeneralRelativity.Tensors.schouten(k.g, k.ricci, k.ricciScalar)));
  sys := Internal.alias(sys, "Schouten", "P", 2);
  sys := Internal.register(sys, "S", 2, Internal.flatten2(GeneralRelativity.Tensors.traceFreeRicci(k.g, k.ricci, k.ricciScalar)));
  sys := Internal.register(sys, "T", 2, Internal.flatten2(
    GeneralRelativity.Tensors.stressEnergyFromEinstein(k.einstein, k.g, gravitationalConstant, c, Lambda)));
  sys := Internal.register(sys, "K", 0, {k.kretschmann});
  annotation(Documentation(info="<html>
<p>The names and conventions are listed in the
<a href=\"modelica://GeneralRelativity.IndexNotation\">package documentation</a>.
The Christoffel symbols are stored as <code>g<sub>il</sub>&Gamma;<sup>l</sup><sub>jk</sub></code>
so that <code>Gamma^i_{jk}</code> in an expression returns
<code>k.Gamma[i, j, k]</code>; they transform as a tensor only in the
expressions where the difference of two connections appears, which is the
user's responsibility.</p>
</html>"));
end systemFromCurvature;
