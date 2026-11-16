within GeneralRelativity.IndexNotation.Internal;
function parseFactor "Read one factor: a real literal, a parenthesized expression or a tensor reference"
  input GeneralRelativity.IndexNotation.TensorSystem sys;
  input String s;
  input Integer start;
  output GeneralRelativity.IndexNotation.TensorValue v;
  output Integer pos;
protected
  String c;
  Integer next "Position returned by a sub-parser; kept apart from pos because some code generators alias an output with a same-named input";
  Real number;
algorithm
  pos := skipSpace(s, start);
  c := charAt(s, pos);
  if c == "(" then
    (v, next) := parseExpression(sys, s, pos + 1);
    pos := skipSpace(s, next);
    if charAt(s, pos) <> ")" then
      fail(s, pos, "expected ')'");
    end if;
    pos := pos + 1;
  elseif isDigit(c) or c == "." then
    (next, number) := Modelica.Utilities.Strings.Advanced.scanReal(s, pos, unsigned = true);
    if next == pos then
      fail(s, pos, "malformed number");
    end if;
    v := scalarValue(number);
    pos := next;
  elseif isLetter(c) or c == "\\" then
    (v, next) := parseTensor(sys, s, pos);
    pos := next;
  elseif c == "" then
    fail(s, pos, "unexpected end of expression");
  else
    fail(s, pos, "unexpected character '" + c + "'");
  end if;
end parseFactor;
