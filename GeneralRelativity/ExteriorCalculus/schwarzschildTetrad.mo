within GeneralRelativity.ExteriorCalculus;
function schwarzschildTetrad "Static orthonormal tetrad of Schwarzschild in {ct, r, theta, phi}; p = {rs}"
  extends TetradField;
protected
  Real rs = p[1] "Schwarzschild radius";
  Real r = x[2];
  Real f = 1.0 - rs / r;
algorithm
  e := zeros(4, 4);
  e[1, 1] := 1.0 / sqrt(f);
  e[2, 2] := sqrt(f);
  e[3, 3] := 1.0 / r;
  e[4, 4] := 1.0 / (r * sin(x[3]));
  annotation(Documentation(info="<html>
<p>Frame of static observers: <code>e<sub>0</sub> = f<sup>-1/2</sup>&part;<sub>0</sub></code>,
<code>e<sub>1</sub> = f<sup>1/2</sup>&part;<sub>r</sub></code>,
<code>e<sub>2</sub> = r<sup>-1</sup>&part;<sub>&theta;</sub></code>,
<code>e<sub>3</sub> = (r sin&theta;)<sup>-1</sup>&part;<sub>&phi;</sub></code>,
<code>f = 1 - r<sub>s</sub>/r</code>, valid outside the horizon. Its coframe is
<code>&theta;<sup>0</sup> = f<sup>1/2</sup>dx<sup>0</sup></code>,
<code>&theta;<sup>1</sup> = f<sup>-1/2</sup>dr</code>, <code>&theta;<sup>2</sup> = r d&theta;</code>,
<code>&theta;<sup>3</sup> = r sin&theta; d&phi;</code> and the nonzero connection forms are
<code>&omega;<sup>0</sup><sub>1</sub> = (f'/2) dx<sup>0</sup></code>,
<code>&omega;<sup>2</sup><sub>1</sub> = f<sup>1/2</sup>d&theta;</code>,
<code>&omega;<sup>3</sup><sub>1</sub> = f<sup>1/2</sup>sin&theta; d&phi;</code>,
<code>&omega;<sup>3</sup><sub>2</sub> = cos&theta; d&phi;</code> and their antisymmetric
partners. Pass with <a href=\"modelica://GeneralRelativity.Metrics.schwarzschild\">Metrics.schwarzschild</a>
and the same <code>p</code>.</p>
</html>"));
end schwarzschildTetrad;
