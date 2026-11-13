within GeneralRelativity.StressEnergy;
function fieldStrength "Covariant Faraday tensor F_ab in an inertial Cartesian frame {ct, x, y, z} from E and B"
  input Real E[3](each unit = "V/m") "Electric field (SI); in Gaussian units pass c = 1 and E in statvolt/cm";
  input Real B[3](each unit = "T") "Magnetic field";
  input Real c(unit = "m/s") = 1.0 "Speed of light; 1 in geometric or Gaussian units";
  output Real F[4, 4] "F_ab, antisymmetric";
algorithm
  F := zeros(4, 4);
  for i in 1:3 loop
    F[1, i + 1] := -E[i] / c;
    F[i + 1, 1] := E[i] / c;
  end for;
  F[2, 3] := B[3];
  F[3, 2] := -B[3];
  F[3, 4] := B[1];
  F[4, 3] := -B[1];
  F[4, 2] := B[2];
  F[2, 4] := -B[2];
  annotation(Documentation(info="<html>
<p>With <code>A<sub>a</sub> = (-&phi;/c, <b>A</b>)</code> and
<code>F<sub>ab</sub> = &part;<sub>a</sub>A<sub>b</sub> - &part;<sub>b</sub>A<sub>a</sub></code>
(<code>x[1] = ct</code>), the components in an inertial Cartesian frame are
<code>F<sub>0i</sub> = -E<sub>i</sub>/c</code> and
<code>F<sub>ij</sub> = &epsilon;<sub>ijk</sub>B<sub>k</sub></code>. The invariants are
<code>F<sub>ab</sub>F<sup>ab</sup> = 2(B&sup2; - E&sup2;/c&sup2;)</code>. Pass the result to
<a href=\"modelica://GeneralRelativity.StressEnergy.electromagnetic\">electromagnetic</a>
with the Minkowski metric.</p>
</html>"));
end fieldStrength;
