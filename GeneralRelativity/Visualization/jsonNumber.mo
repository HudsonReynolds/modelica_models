within GeneralRelativity.Visualization;
function jsonNumber "Real formatted for JSON, 15 significant digits by default"
  input Real value;
  input Integer significantDigits = 15 "Digits kept; large fields use fewer to keep files small";
  output String text;
algorithm
  text := String(value, significantDigits = significantDigits);
end jsonNumber;
