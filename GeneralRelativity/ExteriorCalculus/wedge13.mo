within GeneralRelativity.ExteriorCalculus;
function wedge13 "Wedge product of a 1-form and a 3-form, a 4-form"
  input Real a[4] "1-form";
  input Real H[4, 4, 4] "3-form (totally antisymmetric)";
  output Real V[4, 4, 4, 4] "4-form (a ^ H)_ijkl = a_i H_jkl - a_j H_ikl + a_k H_ijl - a_l H_ijk";
algorithm
  for i in 1:4 loop
    for j in 1:4 loop
      for k in 1:4 loop
        for l in 1:4 loop
          V[i, j, k, l] := a[i] * H[j, k, l] - a[j] * H[i, k, l] + a[k] * H[i, j, l] - a[l] * H[i, j, k];
        end for;
      end for;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p><code>H &and; a = -a &and; H</code>. With <code>H = *b</code>,
<code>a &and; *b = g<sup>-1</sup>(a, b) &epsilon;</code>.</p>
</html>"));
end wedge13;
