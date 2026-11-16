within GeneralRelativity.IndexNotation.Internal;
function expandToken "Index labels of one token: a Greek letter name is one label, anything else one label per character"
  input String token;
  input String s "Whole expression, for error messages";
  input Integer pos "Position of the token, for error messages";
  output String labels[MaxRank];
  output Integer count;
protected
  Integer n = Modelica.Utilities.Strings.length(token);
  String c;
algorithm
  labels := fill("", MaxRank);
  count := 0;
  if isGreek(token) then
    count := 1;
    labels[1] := token;
  else
    for i in 1:n loop
      c := charAt(token, i);
      if not (isLetter(c) or isDigit(c)) then
        fail(s, pos, "index '" + c + "' is not a letter or a digit");
      end if;
      if count >= MaxRank then
        fail(s, pos, "more than " + String(MaxRank) + " indices on one tensor");
      end if;
      count := count + 1;
      labels[count] := c;
    end for;
  end if;
end expandToken;
