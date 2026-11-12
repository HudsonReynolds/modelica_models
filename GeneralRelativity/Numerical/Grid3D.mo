within GeneralRelativity.Numerical;
record Grid3D "Uniform Cartesian grid: node (i, j, k) is at origin + (i-1, j-1, k-1) .* spacing"
  Integer n[3] = {16, 16, 16} "Node count per axis";
  Real origin[3](each unit = "m") = {0.0, 0.0, 0.0} "Position of node (1, 1, 1)";
  Real spacing[3](each unit = "m") = {1.0, 1.0, 1.0} "Node spacing per axis";
end Grid3D;
