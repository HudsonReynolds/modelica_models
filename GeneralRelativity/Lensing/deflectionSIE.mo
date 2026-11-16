within GeneralRelativity.Lensing;
function deflectionSIE "Scaled deflection of a singular isothermal ellipse (major axis along x); q = 1 is the isothermal sphere"
  input Real theta[2](each unit = "rad") "Image-plane position relative to the lens centre, in the lens frame";
  input Real b(unit = "rad") "Deflection scale; the Einstein radius when q = 1";
  input Real q = 1.0 "Axis ratio of the iso-density contours, 0 < q <= 1";
  output Real alpha[2](each unit = "rad");
protected
  Real s = sqrt(max(1.0 - q^2, 0.0));
  Real psi = sqrt(q^2 * theta[1]^2 + theta[2]^2);
algorithm
  if psi <= 0.0 then
    alpha := zeros(2);
  elseif s < 1.0e-6 then
    alpha := b * theta / psi;
  else
    alpha := {b * q / s * atan(s * theta[1] / psi), b * q / s * 0.5 * log((psi + s * theta[2]) / (psi - s * theta[2]))};
  end if;
  annotation(Documentation(info="<html>
<p>Convergence <code>&kappa; = q b/(2&psi;)</code>, <code>&psi; = &radic;(q&sup2;x&sup2; + y&sup2;)</code>
(<a href=\"modelica://GeneralRelativity.Lensing.convergenceSIE\">convergenceSIE</a>), deflection
<code>&alpha;<sub>x</sub> = (bq/s) atan(s x/&psi;)</code>, <code>&alpha;<sub>y</sub> = (bq/s) atanh(s y/&psi;)</code>,
<code>s = &radic;(1 &minus; q&sup2;)</code> (Kormann, Schneider and Bartelmann 1994; Keeton 2001, without core).
The tests check that half its divergence is the convergence and its curl vanishes. Rotate
<code>&theta;</code> into the lens frame for another position angle.</p>
</html>"));
end deflectionSIE;
