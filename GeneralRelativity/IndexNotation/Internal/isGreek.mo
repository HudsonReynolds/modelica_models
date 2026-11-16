within GeneralRelativity.IndexNotation.Internal;
function isGreek "true when a token is the name of a Greek letter and therefore one index, such as mu or nu"
  input String token;
  output Boolean result;
protected
  constant String names[24] = {"alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta",
    "iota", "kappa", "lambda", "mu", "nu", "xi", "omicron", "pi", "rho", "sigma", "tau", "upsilon",
    "phi", "chi", "psi", "omega"};
algorithm
  result := false;
  for i in 1:size(names, 1) loop
    if token == names[i] then
      result := true;
    end if;
  end for;
end isGreek;
