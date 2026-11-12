within GeneralRelativity.Examples;
model SchwarzschildCurvature "Curvature tensors of the Schwarzschild metric at three radii against closed forms"
  parameter Real rs = 2.0 "Schwarzschild radius in geometric units";
  parameter Real step[4] = {1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3} "Finite-difference steps";
  parameter Real radii[3] = {3.0, 5.0, 10.0} * rs "Sample radii";
  final parameter Real theta = GeneralRelativity.Constants.pi / 2.0;
  GeneralRelativity.Tensors.Curvature k1 = GeneralRelativity.Metrics.curvatureAt(
    GeneralRelativity.Metrics.schwarzschild, {0.0, radii[1], theta, 0.0}, {rs}, step);
  GeneralRelativity.Tensors.Curvature k2 = GeneralRelativity.Metrics.curvatureAt(
    GeneralRelativity.Metrics.schwarzschild, {0.0, radii[2], theta, 0.0}, {rs}, step);
  GeneralRelativity.Tensors.Curvature k3 = GeneralRelativity.Metrics.curvatureAt(
    GeneralRelativity.Metrics.schwarzschild, {0.0, radii[3], theta, 0.0}, {rs}, step);
  Real christoffelError[3] "Largest deviation of the numerical Christoffel symbols from the closed form";
  Real ricciNorm[3] "Frobenius norm of the Ricci tensor, zero in vacuum";
  Real kretschmannRelativeError[3] "Relative error against 12 rs^2 / r^6";
  Real impliedStressEnergy[4, 4] "c^4/(8 pi G) G_ab at the innermost radius, zero in vacuum";
  Real tidal[3, 3] "Tidal tensor c^2 R_0i0j at the innermost radius in the coordinate basis";
equation
  christoffelError = {
    maxAbs3(k1.Gamma - GeneralRelativity.Metrics.schwarzschildChristoffel(k1.x, rs)),
    maxAbs3(k2.Gamma - GeneralRelativity.Metrics.schwarzschildChristoffel(k2.x, rs)),
    maxAbs3(k3.Gamma - GeneralRelativity.Metrics.schwarzschildChristoffel(k3.x, rs))};
  ricciNorm = {frobenius(k1.ricci), frobenius(k2.ricci), frobenius(k3.ricci)};
  kretschmannRelativeError = {
    abs(k1.kretschmann / (12.0 * rs^2 / radii[1]^6) - 1.0),
    abs(k2.kretschmann / (12.0 * rs^2 / radii[2]^6) - 1.0),
    abs(k3.kretschmann / (12.0 * rs^2 / radii[3]^6) - 1.0)};
  impliedStressEnergy = GeneralRelativity.Tensors.stressEnergyFromEinstein(k1.einstein, k1.g, 1.0, 1.0);
  tidal = GeneralRelativity.Perturbation.tidalTensor(k1.riemann, 1.0);
  for i in 1:3 loop
    assert(christoffelError[i] < 1.0e-5, "Numerical Christoffel symbols disagree with the closed form");
    assert(ricciNorm[i] < 1.0e-5, "Schwarzschild Ricci tensor is not zero");
    assert(kretschmannRelativeError[i] < 1.0e-4, "Kretschmann scalar disagrees with 12 rs^2 / r^6");
  end for;
  assert(k1.valid > 0.5 and k2.valid > 0.5 and k3.valid > 0.5, "Metric inversion failed");
protected
  function maxAbs3 "Largest absolute entry of a rank-3 array"
    input Real values[:, :, :];
    output Real result;
  algorithm
    result := 0.0;
    for i in 1:size(values, 1) loop
      for j in 1:size(values, 2) loop
        for k in 1:size(values, 3) loop
          result := max(result, abs(values[i, j, k]));
        end for;
      end for;
    end for;
  end maxAbs3;

  function frobenius "Frobenius norm of a matrix"
    input Real values[:, :];
    output Real result;
  algorithm
    result := sqrt(sum(values[i, j]^2 for i in 1:size(values, 1), j in 1:size(values, 2)));
  end frobenius;
  annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 1.0e-8, Interval = 0.1),
    Documentation(info="<html>
<p>The tidal tensor in the coordinate basis at <code>r = 3r<sub>s</sub></code> has
<code>c&sup2;R<sub>0r0r</sub> = -r<sub>s</sub>/r&sup3; &middot; g<sub>rr</sub> &middot; (-g<sub>tt</sub>)</code>
radially and positive transverse entries scaled by <code>r&sup2;</code>: stretching
along the radius, squeezing across it, traceless after raising with the metric.</p>
</html>"));
end SchwarzschildCurvature;
