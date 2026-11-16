within GeneralRelativity.IndexNotation.Internal;
function parseIndexGroup "Read the index labels after a _ or ^: either a braced group of tokens or one bare token"
  input String s;
  input Integer start "Position just after the _ or ^";
  output String labels[MaxRank];
  output Integer count;
  output Integer pos "Position after the group";
protected
  String c;
  String token;
  Integer tokenStart;
  String part[MaxRank];
  Integer partCount;
  Boolean closed;
algorithm
  labels := fill("", MaxRank);
  count := 0;
  pos := start;
  c := charAt(s, pos);
  if c == "{" then
    pos := pos + 1;
    token := "";
    tokenStart := pos;
    closed := false;
    while not closed loop
      c := charAt(s, pos);
      if c == "" then
        fail(s, start, "missing '}'");
      end if;
      if c == "}" or c == " " or c == "," or c == "\\" or c == "\t" or c == "\n" or c == "\r" then
        if token <> "" then
          (part, partCount) := expandToken(token, s, tokenStart);
          for i in 1:partCount loop
            if count >= MaxRank then
              fail(s, tokenStart, "more than " + String(MaxRank) + " indices on one tensor");
            end if;
            count := count + 1;
            labels[count] := part[i];
          end for;
          token := "";
        end if;
        closed := c == "}";
        pos := pos + 1;
        tokenStart := pos;
      else
        token := token + c;
        pos := pos + 1;
      end if;
    end while;
  else
    if c == "\\" then
      pos := pos + 1;
      c := charAt(s, pos);
    end if;
    tokenStart := pos;
    token := "";
    while isLetter(c) or isDigit(c) loop
      token := token + c;
      pos := pos + 1;
      c := charAt(s, pos);
    end while;
    if token == "" then
      fail(s, pos, "expected an index after '" + charAt(s, start - 1) + "'");
    end if;
    (labels, count) := expandToken(token, s, tokenStart);
  end if;
  if count == 0 then
    fail(s, start, "empty index group");
  end if;
end parseIndexGroup;
