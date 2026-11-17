within GeneralRelativity.Cosmology;
function sphericalBessel "Spherical Bessel functions j_0(x) .. j_lMax(x) by Miller's downward recurrence"
  input Integer lMax "Highest order";
  input Real x "Argument, x >= 0";
  output Real j[lMax + 1] "j[l + 1] = j_l(x)";
protected
  Integer L;
  Real jp1;
  Real jl;
  Real jm1;
  Real j0;
  Real j1;
  Real s;
algorithm
  j := zeros(lMax + 1);
  if x < 1.0e-8 then
    j[1] := 1.0;
    return;
  end if;
  // Start well above both l and x: the recurrence downward is stable for every l (the regular solution
  // dominates going down), then normalize with j_0 or j_1, whichever is not near a zero.
  L := integer(max(lMax, x) + 20.0 + sqrt(40.0 * max(lMax, x)));
  jp1 := 0.0;
  jl := 1.0e-300;
  for l in L:-1:1 loop
    if l <= lMax then
      j[l + 1] := jl;
    end if;
    jm1 := (2 * l + 1) / x * jl - jp1;
    jp1 := jl;
    jl := jm1;
    if abs(jl) > 1.0e250 then
      jp1 := jp1 * 1.0e-250;
      jl := jl * 1.0e-250;
      for i in 1:lMax + 1 loop
        j[i] := j[i] * 1.0e-250;
      end for;
    end if;
  end for;
  j[1] := jl;
  j0 := sin(x) / x;
  j1 := sin(x) / x^2 - cos(x) / x;
  s := if abs(sin(x)) > 0.5 then j0 / jl else j1 / jp1;
  j := s * j;
  annotation(Documentation(info="<html>
<p>Recurrence <code>j<sub>l&minus;1</sub> = (2l+1)/x j<sub>l</sub> &minus; j<sub>l+1</sub></code> started at
<code>L &asymp; max(l<sub>max</sub>, x) + &radic;(40 max(l<sub>max</sub>, x))</code> with arbitrary values, rescaled
to avoid overflow, and normalized with <code>j<sub>0</sub> = sin x/x</code> (or <code>j<sub>1</sub></code> near the
zeros of <code>sin x</code>). Accurate to about <code>10<sup>&minus;13</sup></code> for <code>l, x</code> up to several
thousand. These are the radial functions of plane waves, <code>e<sup>ikr&mu;</sup> = &sum;(2l+1)i<sup>l</sup>j<sub>l</sub>(kr)P<sub>l</sub>(&mu;)</code>:
they project a Fourier mode at the last-scattering distance onto multipole <code>l</code> of the sky.</p>
</html>"));
end sphericalBessel;
