within GeneralRelativity.StressEnergy;
function dust "Pressureless non-interacting matter T_ab = rho u_a u_b"
  input Real g[4, 4] "Covariant metric";
  input Real u[4] "Contravariant four-velocity, g(u,u) = -c^2";
  input Real rho(unit = "kg/m3") "Rest-frame mass density";
  output Real T[4, 4] "T_ab";
protected
  Real uLower[4] = g * u;
algorithm
  for a in 1:4 loop
    for b in 1:4 loop
      T[a, b] := rho * uLower[a] * uLower[b];
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>With <code>u = (c, 0, 0, 0)</code> in Minkowski space this gives
<code>T<sub>11</sub> = &rho;c&sup2;</code>, the rest energy density, and nothing
else. Geodesic dust streams are the <i>non-interacting</i> source of
cosmology and of the Oppenheimer-Snyder collapse.</p>
</html>"));
end dust;
