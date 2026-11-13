within GeneralRelativity.Harmonics;
function wigner3j "Wigner 3j symbol (j1 j2 j3; m1 m2 m3); integer or half-integer arguments"
  input Real j1 "First angular momentum";
  input Real j2 "Second angular momentum";
  input Real j3 "Third angular momentum";
  input Real m1 "Projection of j1";
  input Real m2 "Projection of j2";
  input Real m3 "Projection of j3";
  output Real w "(j1 j2 j3; m1 m2 m3); zero when a selection rule fails";
protected
  Integer phase = integer(floor(2.0 * (j1 - j2 - m3) + 0.5)) "2 (j1 - j2 - m3), even for allowed arguments";
algorithm
  w := if mod(phase, 2) <> 0 then 0.0
    else (if mod(div(phase, 2), 2) == 0 then 1.0 else -1.0) / sqrt(2.0 * j3 + 1.0) * clebschGordan(j1, m1, j2, m2, j3, -m3);
  annotation(Documentation(info="<html>
<p>Arguments follow the two-row layout of the symbol: the three angular momenta,
then the three projections. Defined from
<a href=\"modelica://GeneralRelativity.Harmonics.clebschGordan\">clebschGordan</a> by
<code>(j1 j2 j3; m1 m2 m3) = (-1)<sup>j1 - j2 - m3</sup> &lang;j1 m1 j2 m2|j3 -m3&rang;/sqrt(2j3 + 1)</code>.
It vanishes unless <code>m1 + m2 + m3 = 0</code>, is invariant under cyclic permutations of
its columns and changes by <code>(-1)<sup>j1 + j2 + j3</sup></code> under an odd permutation
or under <code>m &rarr; -m</code>.</p>
</html>"));
end wigner3j;
