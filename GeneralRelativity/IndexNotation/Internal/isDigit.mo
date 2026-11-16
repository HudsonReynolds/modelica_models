within GeneralRelativity.IndexNotation.Internal;
function isDigit "true for a single character 0-9"
  input String c;
  output Boolean result;
algorithm
  result := c == "0" or c == "1" or c == "2" or c == "3" or c == "4" or c == "5" or c == "6" or c == "7" or c == "8" or c == "9";
end isDigit;
