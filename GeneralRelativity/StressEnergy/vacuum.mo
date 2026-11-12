within GeneralRelativity.StressEnergy;
function vacuum "T_ab = 0"
  output Real T[4, 4];
algorithm
  T := zeros(4, 4);
  annotation(Inline = true);
end vacuum;
