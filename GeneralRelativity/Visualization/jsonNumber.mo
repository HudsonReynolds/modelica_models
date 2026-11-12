within GeneralRelativity.Visualization;
function jsonNumber "Real formatted for JSON with 15 significant digits"
  input Real value;
  output String text;
algorithm
  text := String(value, significantDigits = 15);
end jsonNumber;
