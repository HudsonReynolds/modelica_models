within GeneralRelativity.ExteriorCalculus;
function hodge0 "Hodge dual of a 0-form: *f = f eps"
  input Real f "0-form";
  input Real g[4, 4] "Covariant metric";
  output Real V[4, 4, 4, 4] "4-form f eps_abcd";
algorithm
  V := f * leviCivitaTensor(g);
  annotation(Documentation(info="<html>
<p><code>*1</code> is the volume form; see
<a href=\"modelica://GeneralRelativity.ExteriorCalculus\">ExteriorCalculus</a> for the
Hodge convention.</p>
</html>"));
end hodge0;
