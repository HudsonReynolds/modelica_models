within GeneralRelativity.Numerical;
function monopoleBoundary "Dirichlet values -GM/|x - centre| on every boundary node of a grid; interior nodes are zero"
  input Grid3D grid;
  input Real GM(unit = "m3/s2") "Gravitational parameter of the enclosed mass";
  input Real centre[3](each unit = "m") "Centre of mass";
  output Real phi[grid.n[1], grid.n[2], grid.n[3]](each unit = "m2/s2");
protected
  Real x[3];
  Real d[3];
algorithm
  phi := zeros(grid.n[1], grid.n[2], grid.n[3]);
  for i in 1:grid.n[1] loop
    for j in 1:grid.n[2] loop
      for k in 1:grid.n[3] loop
        if i == 1 or i == grid.n[1] or j == 1 or j == grid.n[2] or k == 1 or k == grid.n[3] then
          x := gridPoint(grid, i, j, k);
          d := x - centre;
          phi[i, j, k] := -GM / sqrt(d * d);
        end if;
      end for;
    end for;
  end for;
end monopoleBoundary;
