within GeneralRelativity.StressEnergy;
function perfectFluid "Perfect fluid T_ab = (rho + p/c^2) u_a u_b + p g_ab"
  input Real g[4, 4] "Covariant metric";
  input Real u[4] "Contravariant four-velocity, g(u,u) = -c^2";
  input Real rho(unit = "kg/m3") "Rest-frame mass density (total energy density / c^2)";
  input Real pressure(unit = "Pa") "Isotropic pressure";
  input Real c(unit = "m/s") = 1.0 "Speed of light; 1 in geometric units";
  output Real T[4, 4] "T_ab";
protected
  Real uLower[4] = g * u;
algorithm
  for a in 1:4 loop
    for b in 1:4 loop
      T[a, b] := (rho + pressure / c^2) * uLower[a] * uLower[b] + pressure * g[a, b];
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>In the fluid rest frame of Minkowski space the components are
<code>diag(&rho;c&sup2;, p, p, p)</code>. The Tolman-Oppenheimer-Volkoff equations of
<a href=\"modelica://GeneralRelativity.Examples.StarInterior\">Examples.StarInterior</a>
are the static spherically symmetric Einstein equations with this source.</p>
</html>"));
end perfectFluid;
