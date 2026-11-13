within GeneralRelativity.Frames;
function energyConditions "Rest-frame energy density, principal pressures and null/weak/strong/dominant energy conditions of a type I stress-energy tensor"
  input Real g[4, 4] "Covariant metric (g^00 < 0)";
  input Real T[4, 4] "Covariant stress-energy T_ab (x[1] = ct, so T_00 in a rest frame is the energy density)";
  input Real tolerance = 1.0e-10 "Relative tolerance of the inequalities (scaled by the largest of |rho|, |p_i|)";
  output Boolean nullCondition "NEC: rho + p_i >= 0 for all i";
  output Boolean weakCondition "WEC: NEC and rho >= 0";
  output Boolean strongCondition "SEC: NEC and rho + p_1 + p_2 + p_3 >= 0";
  output Boolean dominantCondition "DEC: rho >= |p_i| for all i";
  output Real rho "Energy density in the rest frame (eigenvalue -rho of T^a_b on the timelike eigenvector)";
  output Real p[3] "Principal pressures, ascending";
  output Boolean isTypeI "true when T has a timelike eigenvector (Hawking-Ellis type I) and the iteration converged";
  output Real u[4] "Contravariant unit timelike eigenvector (rest-frame four-velocity, g(u,u) = -1)";
protected
  constant Integer maxSweeps = 60;
  Real E[4, 4] "Current orthonormal frame, rows are vectors";
  Real Tf[4, 4] "Frame components of T";
  Real ei[4];
  Real ej[4];
  Real scale;
  Real off;
  Real t;
  Real phi;
  Real theta;
  Real cs;
  Real sn;
  Real tolAbs;
  Integer sweep;
  Boolean converged;
  Boolean failed;
algorithm
  E := orthonormalTetrad(g);
  converged := false;
  failed := false;
  sweep := 0;
  while not converged and not failed and sweep < maxSweeps loop
    sweep := sweep + 1;
    Tf := E * T * transpose(E);
    scale := max(abs(Tf[i, j]) for i in 1:4, j in 1:4);
    off := 0.0;
    for i in 1:4 loop
      for j in 1:4 loop
        if i <> j then
          off := max(off, abs(Tf[i, j]));
        end if;
      end for;
    end for;
    if off <= 1.0e-14 * scale then
      converged := true;
    else
      for i in 1:3 loop
        for j in i + 1:4 loop
          if not failed then
            Tf := E * T * transpose(E);
            if abs(Tf[i, j]) > 1.0e-15 * scale then
              ei := E[i, :];
              ej := E[j, :];
              if i == 1 then
                // Boost in the (e0, ej) plane: T'_0j = sinh(2 phi)(T_00 + T_jj)/2 + cosh(2 phi) T_0j = 0.
                if abs(2.0 * Tf[1, j]) >= abs(Tf[1, 1] + Tf[j, j]) then
                  failed := true;
                else
                  t := -2.0 * Tf[1, j] / (Tf[1, 1] + Tf[j, j]);
                  phi := 0.25 * log((1.0 + t) / (1.0 - t));
                  E[1, :] := cosh(phi) * ei + sinh(phi) * ej;
                  E[j, :] := sinh(phi) * ei + cosh(phi) * ej;
                end if;
              else
                // Rotation in the (ei, ej) plane: tan(2 theta) = 2 T_ij / (T_ii - T_jj).
                theta := if abs(Tf[i, i] - Tf[j, j]) < 1.0e-300 then 0.25 * GeneralRelativity.Constants.pi * sign(Tf[i, j])
                  else 0.5 * atan(2.0 * Tf[i, j] / (Tf[i, i] - Tf[j, j]));
                cs := cos(theta);
                sn := sin(theta);
                E[i, :] := cs * ei + sn * ej;
                E[j, :] := -sn * ei + cs * ej;
              end if;
            end if;
          end if;
        end for;
      end for;
    end if;
  end while;
  Tf := E * T * transpose(E);
  isTypeI := converged and not failed;
  rho := Tf[1, 1];
  p := {Tf[2, 2], Tf[3, 3], Tf[4, 4]};
  // Sort the three pressures ascending.
  for i in 1:2 loop
    for j in 1:3 - i loop
      if p[j] > p[j + 1] then
        t := p[j];
        p[j] := p[j + 1];
        p[j + 1] := t;
      end if;
    end for;
  end for;
  u := E[1, :];
  tolAbs := tolerance * max(abs(rho), max(abs(p)));
  nullCondition := isTypeI and min(rho + p[i] for i in 1:3) >= -tolAbs;
  weakCondition := nullCondition and rho >= -tolAbs;
  strongCondition := nullCondition and rho + sum(p) >= -tolAbs;
  dominantCondition := isTypeI and rho >= max(abs(p)) - tolAbs;
  annotation(Documentation(info="<html>
<p>A Hawking-Ellis type I tensor has one timelike eigenvector <code>u</code> of
<code>T<sup>a</sup><sub>b</sub></code>, with eigenvalue <code>&minus;&rho;</code>, and three spacelike ones with the
principal pressures <code>p<sub>i</sub></code>; in that frame <code>T<sub>ab</sub> = diag(&rho;, p<sub>1</sub>, p<sub>2</sub>, p<sub>3</sub>)</code>.
All ordinary matter, perfect fluids, radiation and a cosmological constant are type I; null dust and some
fields with null flux are not. The energy conditions then reduce to</p>
<table>
<tr><th>Condition</th><th>Statement</th><th>Type I form</th></tr>
<tr><td>Null (NEC)</td><td><code>T<sub>ab</sub>k<sup>a</sup>k<sup>b</sup> &ge; 0</code> for null <code>k</code></td><td><code>&rho; + p<sub>i</sub> &ge; 0</code></td></tr>
<tr><td>Weak (WEC)</td><td><code>T<sub>ab</sub>v<sup>a</sup>v<sup>b</sup> &ge; 0</code> for timelike <code>v</code></td><td>NEC and <code>&rho; &ge; 0</code></td></tr>
<tr><td>Strong (SEC)</td><td><code>(T<sub>ab</sub> &minus; T g<sub>ab</sub>/2)v<sup>a</sup>v<sup>b</sup> &ge; 0</code></td><td>NEC and <code>&rho; + &Sigma;p<sub>i</sub> &ge; 0</code></td></tr>
<tr><td>Dominant (DEC)</td><td>WEC and <code>&minus;T<sup>a</sup><sub>b</sub>v<sup>b</sup></code> causal</td><td><code>&rho; &ge; |p<sub>i</sub>|</code></td></tr>
</table>
<p><b>Method.</b> Starting from the Eulerian <a href=\"modelica://GeneralRelativity.Frames.orthonormalTetrad\">orthonormal tetrad</a>,
cyclic Jacobi sweeps diagonalize the frame components <code>T<sub>ab</sub></code>: hyperbolic rotations (boosts)
with <code>tanh 2&phi; = &minus;2T<sub>0j</sub>/(T<sub>00</sub> + T<sub>jj</sub>)</code> remove the momentum density and
ordinary rotations remove the shear stresses, until the off-diagonal part is below
<code>10<sup>&minus;14</sup></code> of the largest component. Every step keeps the frame orthonormal, so the
final <code>e<sub>0</sub></code> is the rest frame <code>u</code>.</p>
<p><b>Limits.</b> A boost exists only while <code>|2T<sub>0j</sub>| &lt; |T<sub>00</sub> + T<sub>jj</sub>|</code>; when it fails,
or the sweeps do not converge within 60 iterations, <code>isTypeI</code> is false and all four conditions are
reported false (the other outputs are then the last iterate and not meaningful). Hyperbolic Jacobi iterations are well behaved when
<code>&rho; + p<sub>i</sub></code> has one sign for all <code>i</code> (the pencil <code>T &minus; &lambda;&eta;</code> is definite),
which covers matter satisfying the NEC and its negatives; tensors with <code>&rho; + p<sub>i</sub></code> of mixed sign
may fail to converge and are then reported as not type I. The inequalities are tested with an absolute margin of <code>tolerance</code> times the
largest of <code>|&rho;|</code>, <code>|p<sub>i</sub>|</code>, so a cosmological constant (<code>&rho; + p = 0</code>)
satisfies the NEC marginally.</p>
</html>"));
end energyConditions;
