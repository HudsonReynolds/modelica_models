within GeneralRelativity.Lensing.Inversion;
function syntheticImage "An artificial observation: the rendered image of a known lens plus Gaussian pixel noise from a seeded generator"
  input Real u[:](each unit = "rad");
  input Real v[:](each unit = "rad");
  input Real lens[7] "True lens parameters";
  input Real sources[:, 7] "True sources";
  input Boolean pointMass = false;
  input Real noiseSigma = 0.0 "Standard deviation of the added noise, in brightness units";
  input Integer seed = 1;
  input Integer subpixels = 2;
  input Real psfSigma(unit = "rad") = 0.0;
  output Real image[size(u, 1), size(v, 1)];
protected
  Integer state = mod(abs(seed) * 7919 + 17, 2147483647);
  Real r1;
  Real r2;
algorithm
  image := renderImage(u, v, lens, sources, pointMass, subpixels, psfSigma);
  if noiseSigma > 0.0 then
    for i in 1:size(u, 1) loop
      for j in 1:size(v, 1) loop
        // Park-Miller minimal standard generator (Schrage's product, so 32-bit integers suffice) and a Box-Muller transform.
        state := 16807 * mod(state, 127773) - 2836 * div(state, 127773);
        state := if state <= 0 then state + 2147483647 else state;
        r1 := state / 2147483647.0;
        state := 16807 * mod(state, 127773) - 2836 * div(state, 127773);
        state := if state <= 0 then state + 2147483647 else state;
        r2 := state / 2147483647.0;
        image[i, j] := image[i, j] + noiseSigma * sqrt(-2.0 * log(max(r1, 1.0e-12))) * cos(2.0 * GeneralRelativity.Constants.pi * r2);
      end for;
    end for;
  end if;
end syntheticImage;
