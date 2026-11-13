within GeneralRelativity.Perturbation;
function bransDicke "PPN parameters of Brans-Dicke scalar-tensor theory with coupling omega"
  input Real omega "Brans-Dicke coupling constant; general relativity is the limit omega -> infinity";
  output PPNParameters parameters "gamma = (1 + omega)/(2 + omega), all other parameters at their GR values";
algorithm
  parameters := PPNParameters(gamma = (1.0 + omega) / (2.0 + omega));
  annotation(Documentation(info="<html>
<p>Brans-Dicke theory is fully conservative with <code>&beta; = 1</code>; only
<code>&gamma; = (1 + &omega;)/(2 + &omega;)</code> differs from general relativity, so the
Cassini bound on <code>&gamma;</code> requires <code>&omega; &gt; 4&times;10<sup>4</sup></code>.</p>
</html>"));
end bransDicke;
