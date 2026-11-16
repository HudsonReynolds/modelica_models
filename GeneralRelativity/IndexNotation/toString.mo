within GeneralRelativity.IndexNotation;
function toString "Text form of a value: its free indices and components, for printing with Modelica.Utilities.Streams.print"
  input TensorValue value;
  input Integer significantDigits = 6;
  output String text;
protected
  Integer n = Internal.power4(value.rank);
  Integer idx[Internal.MaxRank];
  String label;
  Boolean any;
algorithm
  text := "[" + Internal.indexList(value) + "]";
  if value.rank == 0 then
    text := text + " = " + String(value.comp[1], significantDigits = significantDigits);
  elseif value.rank == 1 then
    text := text + " = {";
    for i in 1:4 loop
      text := text + (if i > 1 then ", " else "") + String(value.comp[i], significantDigits = significantDigits);
    end for;
    text := text + "}";
  elseif value.rank == 2 then
    text := text + " =\n";
    for i in 1:4 loop
      text := text + (if i == 1 then "  [" else "   ");
      for j in 1:4 loop
        text := text + (if j > 1 then ", " else "") + String(value.comp[(i - 1) * 4 + j], significantDigits = significantDigits);
      end for;
      text := text + (if i < 4 then ";\n" else "]");
    end for;
  else
    // Higher ranks: list the nonzero components with physics labels 0..3.
    any := false;
    for f in 1:n loop
      if value.comp[f] <> 0.0 then
        any := true;
        idx := Internal.decode(f, value.rank);
        label := "";
        for k in 1:value.rank loop
          label := label + (if k > 1 then " " else "") + String(idx[k] - 1);
        end for;
        text := text + "\n  [" + label + "] = " + String(value.comp[f], significantDigits = significantDigits);
      end if;
    end for;
    if not any then
      text := text + " = 0 (all components)";
    end if;
  end if;
  annotation(Documentation(info="<html>
<p>Scalars print as <code>[(scalar)] = value</code>, vectors as a list, matrices
as four rows, and higher ranks as the nonzero components labelled with
coordinate indices <code>0..3</code>. The index list shows each free label with
<code>^</code> for a raised and <code>_</code> for a lowered slot.</p>
</html>"));
end toString;
