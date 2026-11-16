within GeneralRelativity.IndexNotation.Internal;
function charAt "Character i of a string, or the empty string past its end"
  input String s;
  input Integer i;
  output String c;
algorithm
  c := if i >= 1 and i <= Modelica.Utilities.Strings.length(s) then Modelica.Utilities.Strings.substring(s, i, i) else "";
end charAt;
