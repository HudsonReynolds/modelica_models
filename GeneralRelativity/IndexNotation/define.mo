within GeneralRelativity.IndexNotation;
function define "Add or replace a named tensor from flattened components with the given index positions"
  input TensorSystem system;
  input String name "Identifier of letters and digits, starting with a letter";
  input Integer rank "Number of indices, 0..Internal.MaxRank";
  input Real components[:] "At least 4^rank components, row-major with the last index fastest";
  input Boolean up[rank] = fill(false, rank) "true for each contravariant (upper) index of the components as passed";
  output TensorSystem result;
protected
  TensorValue v;
  Integer n = Internal.power4(rank);
  String c;
algorithm
  assert(Modelica.Utilities.Strings.length(name) > 0 and Internal.isLetter(Internal.charAt(name, 1)),
    "IndexNotation.define: the name \"" + name + "\" must start with a letter");
  for i in 2:Modelica.Utilities.Strings.length(name) loop
    c := Internal.charAt(name, i);
    assert(Internal.isLetter(c) or Internal.isDigit(c),
      "IndexNotation.define: the name \"" + name + "\" may contain only letters and digits");
  end for;
  assert(rank >= 0 and rank <= Internal.MaxRank,
    "IndexNotation.define: rank " + String(rank) + " is outside 0.." + String(Internal.MaxRank));
  assert(size(components, 1) >= n,
    "IndexNotation.define: " + name + " needs " + String(n) + " components, got " + String(size(components, 1)));
  v := Internal.emptyValue(rank);
  for i in 1:n loop
    v.comp[i] := components[i];
  end for;
  for k in 1:rank loop
    v.index[k] := String(k);
    v.up[k] := up[k];
    if up[k] then
      v := Internal.setPosition(v, k, false, system.g, system.gInv);
    end if;
  end for;
  result := Internal.register(system, name, rank, v.comp[1:n]);
  annotation(Documentation(info="<html>
<p>The components are lowered with the system's metric before they are
stored, so an expression may use the tensor in any index position. The
convenience wrappers
<a href=\"modelica://GeneralRelativity.IndexNotation.defineScalar\">defineScalar</a>,
<a href=\"modelica://GeneralRelativity.IndexNotation.defineVector\">defineVector</a>,
<a href=\"modelica://GeneralRelativity.IndexNotation.defineMatrix\">defineMatrix</a>,
<a href=\"modelica://GeneralRelativity.IndexNotation.defineTensor3\">defineTensor3</a> and
<a href=\"modelica://GeneralRelativity.IndexNotation.defineTensor4\">defineTensor4</a>
take ordinary arrays. A name may hold one tensor per rank; defining the same
name and rank again replaces the components.</p>
</html>"));
end define;
