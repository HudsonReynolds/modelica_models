within GeneralRelativity.ExteriorCalculus;
function leviCivitaSymbol "Permutation symbol [abcd]: +1 for even permutations of (1,2,3,4), -1 for odd, 0 for repeated indices"
  input Integer a;
  input Integer b;
  input Integer c;
  input Integer d;
  output Integer s "[abcd] with [1234] = +1";
algorithm
  s := sign((b - a) * (c - a) * (d - a) * (c - b) * (d - b) * (d - c));
  annotation(Inline = true, Documentation(info="<html>
<p>The sign of the product of all pairwise differences equals the parity of the
permutation and vanishes for repeated indices. <code>[1234]</code> corresponds
to <code>[0123]</code> in physics index labels.</p>
</html>"));
end leviCivitaSymbol;
