within GeneralRelativity.IndexNotation.Internal;
function parseTensor "Read a tensor name with its index groups, look it up, raise the requested slots, take fixed components and sum repeated labels"
  input GeneralRelativity.IndexNotation.TensorSystem sys;
  input String s;
  input Integer start;
  output GeneralRelativity.IndexNotation.TensorValue v;
  output Integer pos;
protected
  String c;
  String name;
  Integer nameStart;
  String labels[MaxRank];
  Boolean raised[MaxRank];
  Integer count;
  String part[MaxRank];
  Integer partCount;
  Boolean up;
  Integer entry;
  Integer n;
  Integer offset;
  String available;
  Integer repeats;
  Integer k;
  Integer l;
  Boolean found;
  Integer next "Position returned by parseIndexGroup; kept apart from pos because some code generators alias an output with a same-named input";
algorithm
  pos := start;
  if charAt(s, pos) == "\\" then
    pos := pos + 1;
  end if;
  nameStart := pos;
  name := "";
  c := charAt(s, pos);
  while isLetter(c) or isDigit(c) loop
    name := name + c;
    pos := pos + 1;
    c := charAt(s, pos);
  end while;
  if name == "" then
    fail(s, pos, "expected a tensor name");
  end if;
  // Index groups.
  labels := fill("", MaxRank);
  raised := fill(false, MaxRank);
  count := 0;
  while c == "_" or c == "^" loop
    up := c == "^";
    (part, partCount, next) := parseIndexGroup(s, pos + 1);
    pos := next;
    for i in 1:partCount loop
      if count >= MaxRank then
        fail(s, nameStart, "more than " + String(MaxRank) + " indices on " + name);
      end if;
      count := count + 1;
      labels[count] := part[i];
      raised[count] := up;
    end for;
    c := charAt(s, pos);
  end while;
  // Lookup.
  entry := lookup(sys, name, count);
  if entry == 0 then
    available := "";
    for i in 1:sys.count loop
      available := available + (if i > 1 then ", " else "") + sys.name[i] + "(" + String(sys.rank[i]) + ")";
    end for;
    fail(s, nameStart, "no tensor named " + name + " with " + String(count) + " indices; defined (name(indices)): " + available);
  end if;
  v := emptyValue(count);
  n := power4(count);
  offset := sys.offset[entry];
  for i in 1:n loop
    v.comp[i] := sys.data[offset + i - 1];
  end for;
  for i in 1:count loop
    v.index[i] := labels[i];
    if raised[i] then
      v := setPosition(v, i, true, sys.g, sys.gInv);
    end if;
  end for;
  // Fixed numeric components, from the last slot so earlier slot numbers stay valid.
  for i in count:-1:1 loop
    if isDigit(labels[i]) then
      if labels[i] == "0" then
        v := fixSlot(v, i, 1);
      elseif labels[i] == "1" then
        v := fixSlot(v, i, 2);
      elseif labels[i] == "2" then
        v := fixSlot(v, i, 3);
      elseif labels[i] == "3" then
        v := fixSlot(v, i, 4);
      else
        fail(s, nameStart, "index " + labels[i] + " on " + name + " is out of the range 0..3");
      end if;
    end if;
  end for;
  // Repeated labels are summed; three or more copies are an error.
  for i in 1:v.rank loop
    repeats := 0;
    for j in 1:v.rank loop
      if v.index[j] == v.index[i] then
        repeats := repeats + 1;
      end if;
    end for;
    if repeats > 2 then
      fail(s, nameStart, "index " + v.index[i] + " appears " + String(repeats) + " times on " + name);
    end if;
  end for;
  found := true;
  while found loop
    found := false;
    k := 0;
    l := 0;
    for i in 1:v.rank loop
      for j in i + 1:v.rank loop
        if not found and v.index[i] == v.index[j] then
          found := true;
          k := i;
          l := j;
        end if;
      end for;
    end for;
    if found then
      v := contractSlots(v, k, l, sys.g, sys.gInv);
    end if;
  end while;
end parseTensor;
