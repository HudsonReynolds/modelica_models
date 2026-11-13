within GeneralRelativity.Frames;
function weylScalars "Newman-Penrose Weyl scalars Psi0..Psi4 of a Weyl tensor in a null tetrad"
  input Real C[4, 4, 4, 4] "Weyl tensor C_abcd (coordinate components)";
  input Real l[4] "Null vector l";
  input Real n[4] "Null vector n, l.n = -1";
  input Real mRe[4] "Real part of m";
  input Real mIm[4] "Imaginary part of m, m.mbar = 1";
  output Real psiRe[5] "Real parts of Psi0, Psi1, Psi2, Psi3, Psi4";
  output Real psiIm[5] "Imaginary parts of Psi0, Psi1, Psi2, Psi3, Psi4";
protected
  Real F[4, 4, 4, 4] "C on the real basis {l, n, mRe, mIm}";
  // Each NP vector as complex coefficients on {l, n, mRe, mIm}: rows l, n, m, mbar.
  constant Real coefRe[4, 4] = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 1, 0];
  constant Real coefIm[4, 4] = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 1; 0, 0, 0, -1];
  // Slots of Psi_k: C(l,m,l,m), C(l,n,l,m), C(l,m,mbar,n), C(l,n,mbar,n), C(n,mbar,n,mbar).
  constant Integer slots[5, 4] = [1, 3, 1, 3; 1, 2, 1, 3; 1, 3, 4, 2; 1, 2, 4, 2; 2, 4, 2, 4];
  Integer idx[4];
  Real zRe;
  Real zIm;
  Real t;
algorithm
  F := frameComponents4(C, {l, n, mRe, mIm});
  for k in 1:5 loop
    psiRe[k] := 0.0;
    psiIm[k] := 0.0;
    for a in 1:4 loop
      for b in 1:4 loop
        for c in 1:4 loop
          for d in 1:4 loop
            idx := {a, b, c, d};
            zRe := 1.0;
            zIm := 0.0;
            for s in 1:4 loop
              t := zRe * coefRe[slots[k, s], idx[s]] - zIm * coefIm[slots[k, s], idx[s]];
              zIm := zRe * coefIm[slots[k, s], idx[s]] + zIm * coefRe[slots[k, s], idx[s]];
              zRe := t;
            end for;
            psiRe[k] := psiRe[k] + F[a, b, c, d] * zRe;
            psiIm[k] := psiIm[k] + F[a, b, c, d] * zIm;
          end for;
        end for;
      end for;
    end for;
  end for;
  annotation(Documentation(info="<html>
<p>Definitions for the (&minus;+++) signature with <code>l&middot;n = &minus;1</code>, <code>m&middot;m&#772; = 1</code>:</p>
<ul>
<li><code>&Psi;<sub>0</sub> = C<sub>abcd</sub> l<sup>a</sup>m<sup>b</sup>l<sup>c</sup>m<sup>d</sup></code></li>
<li><code>&Psi;<sub>1</sub> = C<sub>abcd</sub> l<sup>a</sup>n<sup>b</sup>l<sup>c</sup>m<sup>d</sup></code></li>
<li><code>&Psi;<sub>2</sub> = C<sub>abcd</sub> l<sup>a</sup>m<sup>b</sup>m&#772;<sup>c</sup>n<sup>d</sup></code></li>
<li><code>&Psi;<sub>3</sub> = C<sub>abcd</sub> l<sup>a</sup>n<sup>b</sup>m&#772;<sup>c</sup>n<sup>d</sup></code></li>
<li><code>&Psi;<sub>4</sub> = C<sub>abcd</sub> n<sup>a</sup>m&#772;<sup>b</sup>n<sup>c</sup>m&#772;<sup>d</sup></code></li>
</ul>
<p>With the package Riemann sign and the radial tetrad of
<a href=\"modelica://GeneralRelativity.Frames.orthonormalTetrad\">orthonormalTetrad</a>/<a href=\"modelica://GeneralRelativity.Frames.nullTetrad\">nullTetrad</a>,
Schwarzschild gives <code>&Psi;<sub>2</sub> = &minus;M/r&sup3; = &minus;r<sub>s</sub>/(2r&sup3;)</code> (geometric units)
and the other four vanish. Far from an isolated source, with <code>l</code> outgoing,
<code>&Psi;<sub>4</sub></code> carries the outgoing gravitational radiation. The scalars depend on the tetrad;
<a href=\"modelica://GeneralRelativity.Frames.weylInvariants\">weylInvariants</a> combines them into tetrad-independent
invariants.</p>
</html>"));
end weylScalars;
