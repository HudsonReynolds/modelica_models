within GeneralRelativity.Symmetry;
function dilation "Dilation (homothety) generator x^a d_a in a Cartesian chart"
  extends GeneralRelativity.Interfaces.VectorFieldFunction;
algorithm
  xi := x;
  annotation(Documentation(info="<html>
<p>In Minkowski space <code>&Lscr;<sub>&xi;</sub>&eta; = 2&eta;</code> and
<code>&nabla;<sub>a</sub>&xi;<sup>a</sup> = 4</code>: a homothetic, hence conformal
Killing, vector that is not Killing. Parameters <code>p</code> are ignored.</p>
</html>"));
end dilation;
