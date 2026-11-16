within GeneralRelativity.IndexNotation.Internal;
function parseExpression "Read a sum of terms with an optional leading sign"
  input GeneralRelativity.IndexNotation.TensorSystem sys;
  input String s;
  input Integer start;
  output GeneralRelativity.IndexNotation.TensorValue v;
  output Integer pos;
protected
  GeneralRelativity.IndexNotation.TensorValue t;
  String c;
  Boolean negative;
  Integer next "Position returned by a sub-parser; kept apart from pos because some code generators alias an output with a same-named input";
algorithm
  pos := skipSpace(s, start);
  c := charAt(s, pos);
  negative := c == "-";
  if c == "-" or c == "+" then
    pos := pos + 1;
  end if;
  (v, next) := parseTerm(sys, s, pos);
  pos := next;
  if negative then
    v := scale(v, -1.0);
  end if;
  pos := skipSpace(s, pos);
  c := charAt(s, pos);
  while c == "+" or c == "-" loop
    (t, next) := parseTerm(sys, s, pos + 1);
    pos := next;
    v := add(v, t, if c == "+" then 1.0 else -1.0, sys.g, sys.gInv);
    pos := skipSpace(s, pos);
    c := charAt(s, pos);
  end while;
end parseExpression;
