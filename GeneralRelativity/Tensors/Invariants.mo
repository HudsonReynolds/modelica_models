within GeneralRelativity.Tensors;
record Invariants "Scalar curvature invariants at one event"
  Real ricciScalar "R";
  Real ricciSquared "R_ab R^ab";
  Real kretschmann "R_abcd R^abcd";
  Real weylSquared "C_abcd C^abcd";
  Real gaussBonnet "Gauss-Bonnet density K - 4 R_ab R^ab + R^2";
  Real cubic[6] "Cubic invariants in the basis of cubicInvariants";
  Real weylCubed "C_ab^cd C_cd^ef C_ef^ab";
  annotation(Documentation(info="<html>
<p>Filled from a <a href=\"modelica://GeneralRelativity.Tensors.Curvature\">Curvature</a> record by
<a href=\"modelica://GeneralRelativity.Tensors.curvatureInvariants\">curvatureInvariants</a>.
In four dimensions <code>weylSquared = kretschmann &minus; 2 ricciSquared + ricciScalar&sup2;/3</code>.
The entries of <code>cubic</code> are listed in
<a href=\"modelica://GeneralRelativity.Tensors.cubicInvariants\">cubicInvariants</a>.</p>
</html>"));
end Invariants;
