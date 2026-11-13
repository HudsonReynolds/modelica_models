within GeneralRelativity.ExteriorCalculus;
function wedge12 "Wedge product of a 1-form and a 2-form, (a ^ F)_ijk = a_i F_jk + a_j F_ki + a_k F_ij"
  input Real a[4] "1-form";
  input Real F[4, 4] "2-form (antisymmetric)";
  output Real H[4, 4, 4] "3-form a ^ F";
algorithm
  for i in 1:4 loop
    for j in 1:4 loop
      for k in 1:4 loop
        H[i, j, k] := a[i] * F[j, k] + a[j] * F[k, i] + a[k] * F[i, j];
      end for;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>For 2-forms <code>F &and; a = a &and; F</code>, so this also gives the
(2,1) product.</p>
</html>"));
end wedge12;
