within GeneralRelativity.IndexNotation.Internal;
function parseTerm "Read a product: factors joined by juxtaposition or *, and divided by scalar factors with /"
  input GeneralRelativity.IndexNotation.TensorSystem sys;
  input String s;
  input Integer start;
  output GeneralRelativity.IndexNotation.TensorValue v;
  output Integer pos;
protected
  GeneralRelativity.IndexNotation.TensorValue t;
  String c;
  Integer at;
  Boolean more;
  Integer next "Position returned by parseFactor; kept apart from pos because some code generators alias an output with a same-named input";
algorithm
  (v, pos) := parseFactor(sys, s, start);
  more := true;
  while more loop
    pos := skipSpace(s, pos);
    c := charAt(s, pos);
    if c == "*" then
      (t, next) := parseFactor(sys, s, pos + 1);
      pos := next;
      v := multiply(v, t, sys.g, sys.gInv);
    elseif c == "/" then
      at := pos;
      (t, next) := parseFactor(sys, s, pos + 1);
      pos := next;
      if t.rank <> 0 then
        fail(s, at, "the divisor has free indices [" + indexList(t) + "]");
      end if;
      if t.comp[1] == 0.0 then
        fail(s, at, "division by zero");
      end if;
      v := scale(v, 1.0 / t.comp[1]);
    elseif isLetter(c) or isDigit(c) or c == "." or c == "(" or c == "\\" then
      (t, next) := parseFactor(sys, s, pos);
      pos := next;
      v := multiply(v, t, sys.g, sys.gInv);
    else
      more := false;
    end if;
  end while;
end parseTerm;
