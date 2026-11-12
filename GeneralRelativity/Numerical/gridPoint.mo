within GeneralRelativity.Numerical;
function gridPoint "Position of node (i, j, k)"
  input Grid3D grid;
  input Integer i;
  input Integer j;
  input Integer k;
  output Real x[3](each unit = "m");
algorithm
  x := grid.origin + {(i - 1) * grid.spacing[1], (j - 1) * grid.spacing[2], (k - 1) * grid.spacing[3]};
  annotation(Inline = true);
end gridPoint;
