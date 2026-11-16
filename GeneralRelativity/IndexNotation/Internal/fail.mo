within GeneralRelativity.IndexNotation.Internal;
function fail "Stop with a message that quotes the expression and the character position"
  input String s "Expression";
  input Integer pos "Character position, 1-based";
  input String message;
algorithm
  assert(false, "IndexNotation: " + message + " at character " + String(pos) + " of \"" + s + "\"");
end fail;
