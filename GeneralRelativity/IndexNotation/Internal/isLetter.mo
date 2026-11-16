within GeneralRelativity.IndexNotation.Internal;
function isLetter "true for a single ASCII letter a-z or A-Z"
  input String c;
  output Boolean result;
protected
  import Modelica.Utilities.Strings.compare;
  import Modelica.Utilities.Types.Compare;
algorithm
  result := Modelica.Utilities.Strings.length(c) == 1 and (
    (compare(c, "a") <> Compare.Less and compare(c, "z") <> Compare.Greater) or
    (compare(c, "A") <> Compare.Less and compare(c, "Z") <> Compare.Greater));
end isLetter;
