within GeneralRelativity.Harmonics;
function logFactorial "Natural logarithm of n! for a non-negative integer n"
  input Integer n "Argument, n >= 0";
  output Real y "log(n!)";
algorithm
  assert(n >= 0, "logFactorial needs a non-negative argument");
  y := 0.0;
  for k in 2:n loop
    y := y + log(k);
  end for;
  annotation(Documentation(info="<html>
<p>A plain sum of logarithms: exact to rounding and fast for the arguments of the
coupling coefficients (<code>n</code> up to a few hundred), and free of the overflow of
<code>n!</code> beyond <code>n = 170</code>.</p>
</html>"));
end logFactorial;
