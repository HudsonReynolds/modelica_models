within GeneralRelativity.Frames;
function classifyPetrovMatrix "Petrov type from the eigenstructure of the complex symmetric traceless matrix Q = E + i B"
  input Real QRe[3, 3] "Re Q, e.g. the electric part E of the Weyl tensor";
  input Real QIm[3, 3] "Im Q, e.g. the magnetic part B of the Weyl tensor";
  input Real scale "Reference magnitude; |Q| <= tolerance * scale is type O";
  input Real tolerance = 1.0e-3 "Relative tolerance of the degeneracy tests";
  output PetrovType petrov "Petrov type";
protected
  Real norm;
  Real A[3, 3] "Re Q / |Q|";
  Real B[3, 3] "Im Q / |Q|";
  Real A2[3, 3] "Re Q^2";
  Real B2[3, 3] "Im Q^2";
  Real aRe "Re of a = tr(Q^2)/2";
  Real aIm;
  Real bRe "Re of b = tr(Q^3)/3 = det Q";
  Real bIm;
  Real dRe "Re of 4 a^3 - 27 b^2";
  Real dIm;
  Real a2Re;
  Real a2Im;
  Real lRe "Re of the repeated eigenvalue -3b/(2a)";
  Real lIm;
  Real MRe[3, 3] "Re of (Q - lambda)(Q + 2 lambda) = Q^2 + lambda Q - 2 lambda^2";
  Real MIm[3, 3];
algorithm
  norm := sqrt(sum(QRe[i, j]^2 + QIm[i, j]^2 for i in 1:3, j in 1:3));
  if norm <= tolerance * scale or norm <= 0.0 then
    petrov := PetrovType.O;
  else
    A := QRe / norm;
    B := QIm / norm;
    A2 := A * A - B * B;
    B2 := A * B + B * A;
    aRe := 0.5 * sum(A2[i, i] for i in 1:3);
    aIm := 0.5 * sum(B2[i, i] for i in 1:3);
    // tr(Q^3) = tr(Q^2 Q).
    bRe := sum(A2[i, j] * A[j, i] - B2[i, j] * B[j, i] for i in 1:3, j in 1:3) / 3.0;
    bIm := sum(A2[i, j] * B[j, i] + B2[i, j] * A[j, i] for i in 1:3, j in 1:3) / 3.0;
    if sqrt(aRe^2 + aIm^2) < tolerance and sqrt(bRe^2 + bIm^2) < tolerance then
      // All eigenvalues zero: Q is nilpotent.
      petrov := if sqrt(sum(A2[i, j]^2 + B2[i, j]^2 for i in 1:3, j in 1:3)) < tolerance then PetrovType.N else PetrovType.III;
    else
      a2Re := aRe^2 - aIm^2;
      a2Im := 2.0 * aRe * aIm;
      dRe := 4.0 * (a2Re * aRe - a2Im * aIm) - 27.0 * (bRe^2 - bIm^2);
      dIm := 4.0 * (a2Re * aIm + a2Im * aRe) - 27.0 * 2.0 * bRe * bIm;
      if sqrt(dRe^2 + dIm^2) > tolerance then
        petrov := PetrovType.I;
      else
        // lambda = -3 b / (2 a) = -1.5 b conj(a) / |a|^2.
        lRe := -1.5 * (bRe * aRe + bIm * aIm) / (aRe^2 + aIm^2);
        lIm := -1.5 * (bIm * aRe - bRe * aIm) / (aRe^2 + aIm^2);
        MRe := A2 + lRe * A - lIm * B - 2.0 * (lRe^2 - lIm^2) * identity(3);
        MIm := B2 + lRe * B + lIm * A - 4.0 * lRe * lIm * identity(3);
        petrov := if sqrt(sum(MRe[i, j]^2 + MIm[i, j]^2 for i in 1:3, j in 1:3)) < tolerance then PetrovType.D else PetrovType.II;
      end if;
    end if;
  end if;
  annotation(Documentation(info="<html>
<p><code>Q</code> is normalized to unit Frobenius norm. Its eigenvalues are the roots of
<code>&lambda;&sup3; &minus; a&lambda; &minus; b = 0</code> with <code>a = tr(Q&sup2;)/2</code> and
<code>b = det Q = tr(Q&sup3;)/3</code> (Cayley-Hamilton for a traceless matrix). Decision tree:</p>
<ol>
<li><code>|Q| &le; tolerance&middot;scale</code>: <b>O</b>;</li>
<li><code>a &asymp; 0</code> and <code>b &asymp; 0</code> (nilpotent): <b>N</b> if <code>Q&sup2; &asymp; 0</code>, else <b>III</b>;</li>
<li>discriminant <code>4a&sup3; &minus; 27b&sup2;</code> not small: three distinct eigenvalues, <b>I</b>;</li>
<li>otherwise a repeated eigenvalue <code>&lambda; = &minus;3b/(2a)</code> (the other is <code>&minus;2&lambda;</code>):
<b>D</b> if <code>(Q &minus; &lambda;)(Q + 2&lambda;) &asymp; 0</code> (diagonalizable), else <b>II</b>.</li>
</ol>
<p>Every &ldquo;&asymp; 0&rdquo; test compares a quantity of the normalized <code>Q</code> with
<code>tolerance</code>. Curvature from finite differences carries relative errors around
<code>10<sup>&minus;5</sup></code>, so the default <code>10<sup>&minus;3</sup></code> separates the types reliably
unless the field is itself close to a special type (then the answer is the nearby special type).</p>
</html>"));
end classifyPetrovMatrix;
