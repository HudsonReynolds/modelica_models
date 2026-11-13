within GeneralRelativity.NewtonCartan;
function compatibilityAt "Largest components of nabla tau and nabla h for a connection: zero for a Newton-Cartan connection"
  input GeneralRelativity.Interfaces.ConnectionFunction connection "Connection to test";
  input Real x[4] "Coordinates";
  input Real p[:] "Connection parameters";
  input Real step[4] "Finite-difference step for d h";
  input Boolean spherical = false "Chart of x: true for {ct, r, theta, phi}";
  output Real clockResidual "max |nabla_a tau_b|";
  output Real cometricResidual "max |nabla_a h^bc|";
protected
  Real Gamma[4, 4, 4];
  Real h[4, 4];
  Real hPlus[4, 4];
  Real hMinus[4, 4];
  Real dh[4, 4, 4];
  Real tau[4] = clockForm();
  Real xs[4];
  Real value;
algorithm
  Gamma := connection(x, p);
  h := spatialCometric(x, spherical);
  for a in 1:4 loop
    xs := x;
    xs[a] := x[a] + step[a];
    hPlus := spatialCometric(xs, spherical);
    xs[a] := x[a] - step[a];
    hMinus := spatialCometric(xs, spherical);
    dh[a, :, :] := (hPlus - hMinus) / (2.0 * step[a]);
  end for;
  clockResidual := 0.0;
  cometricResidual := 0.0;
  for a in 1:4 loop
    for b in 1:4 loop
      // tau is constant, so nabla_a tau_b = -Gamma^c_ab tau_c.
      clockResidual := max(clockResidual, abs(sum(Gamma[c, a, b] * tau[c] for c in 1:4)));
      for c in 1:4 loop
        value := dh[a, b, c] + sum(Gamma[b, a, d] * h[d, c] + Gamma[c, a, d] * h[b, d] for d in 1:4);
        cometricResidual := max(cometricResidual, abs(value));
      end for;
    end for;
  end for;
end compatibilityAt;
