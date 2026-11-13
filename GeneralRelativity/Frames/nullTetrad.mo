within GeneralRelativity.Frames;
function nullTetrad "Newman-Penrose null tetrad l, n, m = mRe + i mIm from an orthonormal tetrad"
  input Real e[4, 4] "Orthonormal tetrad, e[a, :] = e_a^mu, e[1, :] timelike";
  output Real l[4] "Outgoing null vector (e0 + e1)/sqrt(2)";
  output Real n[4] "Ingoing null vector (e0 - e1)/sqrt(2)";
  output Real mRe[4] "Real part of m = (e2 + i e3)/sqrt(2)";
  output Real mIm[4] "Imaginary part of m";
algorithm
  l := (e[1, :] + e[2, :]) / sqrt(2.0);
  n := (e[1, :] - e[2, :]) / sqrt(2.0);
  mRe := e[3, :] / sqrt(2.0);
  mIm := e[4, :] / sqrt(2.0);
  annotation(Documentation(info="<html>
<p>Normalization for the (&minus;+++) signature: <code>l&middot;n = &minus;1</code>,
<code>m&middot;m&#772; = 1</code>, all other products zero (in particular <code>m&middot;m = 0</code>,
i.e. <code>|mRe|&sup2; = |mIm|&sup2;</code> and <code>mRe&middot;mIm = 0</code>). With the
<a href=\"modelica://GeneralRelativity.Frames.orthonormalTetrad\">orthonormalTetrad</a> of a spherical chart,
<code>l</code> and <code>n</code> point radially outward and inward.</p>
</html>"));
end nullTetrad;
