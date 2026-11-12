within GeneralRelativity.Perturbation.Potentials;
function circularBodyPosition "Position of a body on a circular orbit in the xy plane: d (cos(Omega t + phi0), sin(Omega t + phi0), 0)"
  input Real t(unit = "s") "Coordinate time";
  input Real d(unit = "m") "Orbit radius";
  input Real Omega(unit = "rad/s") "Orbital rate";
  input Real phi0(unit = "rad") = 0.0 "Phase at t = 0";
  output Real position[3](each unit = "m");
algorithm
  position := d * {cos(Omega * t + phi0), sin(Omega * t + phi0), 0.0};
  annotation(Inline = true);
end circularBodyPosition;
