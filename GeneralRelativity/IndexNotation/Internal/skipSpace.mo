within GeneralRelativity.IndexNotation.Internal;
function skipSpace "Position of the first character at or after i that is not whitespace"
  input String s;
  input Integer i;
  output Integer pos;
protected
  String c;
algorithm
  pos := i;
  c := charAt(s, pos);
  while c == " " or c == "\t" or c == "\n" or c == "\r" loop
    pos := pos + 1;
    c := charAt(s, pos);
  end while;
end skipSpace;
