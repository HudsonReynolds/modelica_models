within GeneralRelativity.Numerical;
function uniformSphereDensity "Mass density of a homogeneous sphere sampled on a grid"
  input Grid3D grid;
  input Real rho(unit = "kg/m3") "Density inside the sphere";
  input Real radius(unit = "m") "Sphere radius";
  input Real centre[3](each unit = "m") "Sphere centre";
  output Real density[grid.n[1], grid.n[2], grid.n[3]](each unit = "kg/m3");
protected
  Real d[3];
algorithm
  for i in 1:grid.n[1] loop
    for j in 1:grid.n[2] loop
      for k in 1:grid.n[3] loop
        d := gridPoint(grid, i, j, k) - centre;
        density[i, j, k] := if d * d <= radius^2 then rho else 0.0;
      end for;
    end for;
  end for;
end uniformSphereDensity;
