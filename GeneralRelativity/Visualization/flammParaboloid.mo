within GeneralRelativity.Visualization;
function flammParaboloid "Height of Flamm's paraboloid z = 2 sqrt(rs (r - rs)) embedding the Schwarzschild equatorial plane"
  input Real r "Areal radius, r >= rs";
  input Real rs "Schwarzschild radius";
  output Real z;
algorithm
  z := 2.0 * sqrt(max(rs * (r - rs), 0.0));
  annotation(Inline = true);
end flammParaboloid;
