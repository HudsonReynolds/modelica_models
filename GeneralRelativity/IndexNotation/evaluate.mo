within GeneralRelativity.IndexNotation;
function evaluate "Evaluate an index expression against a system; the result keeps the free indices"
  input TensorSystem system;
  input String expression "For example \"g_{mu nu} g^{mu nu}\" or \"R_{ab} - 1/2 g_{ab} R\"";
  output TensorValue value;
protected
  Integer pos;
algorithm
  (value, pos) := Internal.parseExpression(system, expression, 1);
  pos := Internal.skipSpace(expression, pos);
  if pos <= Modelica.Utilities.Strings.length(expression) then
    Internal.fail(expression, pos, "unexpected character '" + Internal.charAt(expression, pos) + "'");
  end if;
  annotation(Documentation(info="<html>
<p>See the <a href=\"modelica://GeneralRelativity.IndexNotation\">package documentation</a>
for the syntax. The free indices of the result are ordered by first
appearance in the expression and carry the positions the expression gave
them; <a href=\"modelica://GeneralRelativity.IndexNotation.scalar\">scalar</a>,
<a href=\"modelica://GeneralRelativity.IndexNotation.vector\">vector</a>,
<a href=\"modelica://GeneralRelativity.IndexNotation.matrix\">matrix</a>,
<a href=\"modelica://GeneralRelativity.IndexNotation.tensor3\">tensor3</a> and
<a href=\"modelica://GeneralRelativity.IndexNotation.tensor4\">tensor4</a>
unpack them into arrays.</p>
</html>"));
end evaluate;
