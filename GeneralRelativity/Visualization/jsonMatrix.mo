within GeneralRelativity.Visualization;
function jsonMatrix "Real matrix formatted as a JSON array of row arrays"
  input Real values[:, :];
  input Integer significantDigits = 15 "Digits kept per element";
  output String text;
algorithm
  text := "[";
  for i in 1:size(values, 1) loop
    text := text + (if i > 1 then "," else "") + jsonVector(values[i, :], significantDigits);
  end for;
  text := text + "]";
end jsonMatrix;
