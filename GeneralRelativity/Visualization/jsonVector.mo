within GeneralRelativity.Visualization;
function jsonVector "Real vector formatted as a JSON array"
  input Real values[:];
  input Integer significantDigits = 15 "Digits kept per element";
  output String text;
algorithm
  text := "[";
  for i in 1:size(values, 1) loop
    text := text + (if i > 1 then "," else "") + jsonNumber(values[i], significantDigits);
  end for;
  text := text + "]";
end jsonVector;
