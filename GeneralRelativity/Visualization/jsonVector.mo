within GeneralRelativity.Visualization;
function jsonVector "Real vector formatted as a JSON array"
  input Real values[:];
  output String text;
algorithm
  text := "[";
  for i in 1:size(values, 1) loop
    text := text + (if i > 1 then "," else "") + jsonNumber(values[i]);
  end for;
  text := text + "]";
end jsonVector;
