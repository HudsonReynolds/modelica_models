within GeneralRelativity.Cosmology;
function cmbSpectrum "CMB temperature power spectrum D_l = l(l+1)C_l/(2 pi) in uK^2 by the line-of-sight integral"
  input Parameters par;
  input Real ells[:] "Multipoles (increasing)";
  input Real k[:] "Wavenumbers of the sources (increasing)";
  input Real etaRec[:] "Conformal times of the recombination sources";
  input Real rec[size(k, 1), size(etaRec, 1), 5] "delta_gamma/4, psi, theta_b/k, psi' + phi', Pi (acousticSources)";
  input Real etaIsw[:] "Conformal times of the ISW source";
  input Real isw[size(k, 1), size(etaIsw, 1)] "psi' + phi'";
  input Real gRec[size(etaRec, 1)] "Visibility at etaRec";
  input Real tauRec[size(etaRec, 1)] "Optical depth at etaRec";
  input Real tauIsw[size(etaIsw, 1)] "Optical depth at etaIsw";
  input Real eta0 "Conformal time today";
  input Real kMax = 0.36 "Largest wavenumber of the l integrals";
  input Real dkFine = 4.0e-5 "Fine-grid spacing above kLog";
  input Real kLog = 2.0e-3 "Below this the fine grid is logarithmic";
  input Real dx = 0.4 "Spacing of the Bessel table in x = k (eta0 - eta)";
  input Integer lIswMax = 400 "The ISW term is evaluated up to this l";
  output Real Dl[size(ells, 1)] "Total";
  output Real Dsw[size(ells, 1)] "Monopole (Sachs-Wolfe and acoustic density) term alone";
  output Real Ddop[size(ells, 1)] "Doppler term alone";
  output Real Disw[size(ells, 1)] "Integrated Sachs-Wolfe term alone";
protected
  Integer nl = size(ells, 1);
  Integer nk = size(k, 1);
  Integer nr = size(etaRec, 1);
  Integer ni = size(etaIsw, 1);
  Integer lMax = integer(max(ells));
  Integer nx = integer(kMax * (eta0 - etaRec[1]) / dx) + 3;
  Real J[nl, nx] "j_l at x = (i - 1) dx";
  Real Jd[nl, nx] "j_l'";
  Real jAll[lMax + 1];
  Real x;
  Real H0, Og, On, Ob, Oc, OL;
  Real wr[nr];
  Real wi[ni];
  Real M0[nk, nr] "Spline second derivatives (in k) of the monopole source";
  Real Mp[nk, nr];
  Real M1[nk, nr];
  Real M2[nk, nr];
  Real M3[nk, nr];
  Real f3[nr];
  Real S3[nr] "Polarization (photon quadrupole) source";
  Real jddv;
  Real f2[nr];
  Real S2[nr] "ISW source during recombination";
  Real MI[nk, ni];
  Real f0[nr];
  Real fp[nr];
  Real f1[nr];
  Real fI[ni];
  Real S0[nr];
  Real S1[nr];
  Real SI[ni];
  Integer nLog = 200;
  Integer nLin = integer((kMax - kLog) / dkFine);
  Integer nf = nLog + nLin;
  Real kf;
  Real kPrev;
  Real kNext;
  Real wk;
  Integer jk;
  Real t;
  Real hk;
  Real A;
  Real B;
  Real D;
  Real th0;
  Real th1;
  Real thI;
  Real P;
  Real C0[nl];
  Real Csw[nl];
  Real Cd[nl];
  Real Ci[nl];
  Real l;
  Real xl;
  Real u;
  Integer ix;
  Real tt;
  Real h00, h10, h01, h11;
  Real jv;
  Real jdv;
  Real jddLo;
  Real jddHi;
  Real xLo;
  Real xHi;
  Real xMinL;
algorithm
  (H0, Og, On, Ob, Oc, OL) := densities(par);
  // Bessel table: j_l and j_l' for the requested l on a uniform grid in x (Miller recurrence per node).
  for i in 1:nx loop
    x := max((i - 1) * dx, 1.0e-6);
    jAll := sphericalBessel(lMax, x);
    for n in 1:nl loop
      J[n, i] := jAll[integer(ells[n]) + 1];
      Jd[n, i] := (if ells[n] > 0.5 then jAll[integer(ells[n])] else -jAll[2]) - (if ells[n] > 0.5 then (ells[n] + 1.0) / x * jAll[integer(ells[n]) + 1] else 0.0);
    end for;
  end for;
  // Quadrature weights in eta (trapezoid) and the source factors that do not depend on k.
  for j in 1:nr loop
    wr[j] := 0.5 * (etaRec[min(j + 1, nr)] - etaRec[max(j - 1, 1)]);
  end for;
  for j in 1:ni loop
    wi[j] := 0.5 * (etaIsw[min(j + 1, ni)] - etaIsw[max(j - 1, 1)]) * exp(-tauIsw[j]);
  end for;
  // Natural cubic splines in k of each source column.
  M0 := splineSecondDerivatives(k, rec[:, :, 1]);
  Mp := splineSecondDerivatives(k, rec[:, :, 2]);
  M1 := splineSecondDerivatives(k, rec[:, :, 3]);
  M2 := splineSecondDerivatives(k, rec[:, :, 4]);
  M3 := splineSecondDerivatives(k, rec[:, :, 5]);
  MI := splineSecondDerivatives(k, isw);
  C0 := zeros(nl);
  Csw := zeros(nl);
  Cd := zeros(nl);
  Ci := zeros(nl);
  for q in 1:nf loop
    // Fine wavenumbers: logarithmic up to kLog, then uniform; trapezoid weight wk.
    kf := if q <= nLog then k[1] * (kLog / k[1])^((q - 1) / nLog) else kLog + (q - nLog) * dkFine;
    kPrev := if q == 1 then kf else (if q - 1 <= nLog then k[1] * (kLog / k[1])^((q - 2) / nLog) else kLog + (q - 1 - nLog) * dkFine);
    kNext := if q == nf then kf else (if q + 1 <= nLog then k[1] * (kLog / k[1])^(q / nLog) else kLog + (q + 1 - nLog) * dkFine);
    wk := 0.5 * (kNext - kPrev);
    // Spline evaluation at kf.
    jk := 1;
    while jk < nk - 1 and k[jk + 1] < kf loop
      jk := jk + 1;
    end while;
    hk := k[jk + 1] - k[jk];
    A := (k[jk + 1] - kf) / hk;
    B := 1.0 - A;
    for j in 1:nr loop
      f0[j] := A * rec[jk, j, 1] + B * rec[jk + 1, j, 1] + ((A^3 - A) * M0[jk, j] + (B^3 - B) * M0[jk + 1, j]) * hk^2 / 6.0;
      fp[j] := A * rec[jk, j, 2] + B * rec[jk + 1, j, 2] + ((A^3 - A) * Mp[jk, j] + (B^3 - B) * Mp[jk + 1, j]) * hk^2 / 6.0;
      f1[j] := A * rec[jk, j, 3] + B * rec[jk + 1, j, 3] + ((A^3 - A) * M1[jk, j] + (B^3 - B) * M1[jk + 1, j]) * hk^2 / 6.0;
      S0[j] := gRec[j] * wr[j] * (f0[j] + fp[j]);
      S1[j] := gRec[j] * wr[j] * f1[j];
      f3[j] := A * rec[jk, j, 5] + B * rec[jk + 1, j, 5] + ((A^3 - A) * M3[jk, j] + (B^3 - B) * M3[jk + 1, j]) * hk^2 / 6.0;
      S3[j] := gRec[j] * wr[j] * f3[j] / 16.0;
      f2[j] := A * rec[jk, j, 4] + B * rec[jk + 1, j, 4] + ((A^3 - A) * M2[jk, j] + (B^3 - B) * M2[jk + 1, j]) * hk^2 / 6.0;
      S2[j] := exp(-tauRec[j]) * wr[j] * f2[j];
    end for;
    for j in 1:ni loop
      fI[j] := A * isw[jk, j] + B * isw[jk + 1, j] + ((A^3 - A) * MI[jk, j] + (B^3 - B) * MI[jk + 1, j]) * hk^2 / 6.0;
      SI[j] := wi[j] * fI[j];
    end for;
    P := par.AS * (kf / par.kPivot)^(par.nS - 1.0);
    for n in 1:nl loop
      l := ells[n];
      // j_l(x) is negligible below its turning point region; skip modes that never reach it.
      xMinL := max(0.0, l - 8.0 * l^(1.0 / 3.0) - 5.0);
      th0 := 0.0;
      th1 := 0.0;
      thI := 0.0;
      if kf * (eta0 - etaRec[1]) >= xMinL then
        for j in 1:nr loop
          xl := kf * (eta0 - etaRec[j]);
          if xl >= xMinL then
            u := xl / dx;
            ix := min(integer(u), nx - 2);
            tt := u - ix;
            h00 := (1.0 + 2.0 * tt) * (1.0 - tt)^2;
            h10 := tt * (1.0 - tt)^2;
            h01 := tt^2 * (3.0 - 2.0 * tt);
            h11 := tt^2 * (tt - 1.0);
            xLo := max(ix * dx, 1.0e-6);
            xHi := (ix + 1) * dx;
            jv := h00 * J[n, ix + 1] + h10 * dx * Jd[n, ix + 1] + h01 * J[n, ix + 2] + h11 * dx * Jd[n, ix + 2];
            // j_l'' from Bessel's equation, for the Hermite interpolation of j_l'.
            jddLo := -2.0 / xLo * Jd[n, ix + 1] - (1.0 - l * (l + 1.0) / xLo^2) * J[n, ix + 1];
            jddHi := -2.0 / xHi * Jd[n, ix + 2] - (1.0 - l * (l + 1.0) / xHi^2) * J[n, ix + 2];
            jdv := h00 * Jd[n, ix + 1] + h10 * dx * jddLo + h01 * Jd[n, ix + 2] + h11 * dx * jddHi;
            // j_l'' at x from Bessel's equation, for the polarization term g Pi (j_l + 3 j_l'') / 16.
            jddv := -2.0 / xl * jdv - (1.0 - l * (l + 1.0) / xl^2) * jv;
            th0 := th0 + S0[j] * jv + S3[j] * (jv + 3.0 * jddv);
            th1 := th1 + S1[j] * jdv;
            thI := thI + S2[j] * jv;
          end if;
        end for;
      end if;
      if l <= lIswMax then
        for j in 1:ni loop
          xl := kf * (eta0 - etaIsw[j]);
          if xl >= xMinL then
            u := xl / dx;
            ix := min(integer(u), nx - 2);
            tt := u - ix;
            h00 := (1.0 + 2.0 * tt) * (1.0 - tt)^2;
            h10 := tt * (1.0 - tt)^2;
            h01 := tt^2 * (3.0 - 2.0 * tt);
            h11 := tt^2 * (tt - 1.0);
            thI := thI + SI[j] * (h00 * J[n, ix + 1] + h10 * dx * Jd[n, ix + 1] + h01 * J[n, ix + 2] + h11 * dx * Jd[n, ix + 2]);
          end if;
        end for;
      end if;
      // C_l = 4 pi int dk/k Delta_R^2(k) |Theta_l(k)|^2
      C0[n] := C0[n] + wk / kf * P * (th0 + th1 + thI)^2;
      Csw[n] := Csw[n] + wk / kf * P * th0^2;
      Cd[n] := Cd[n] + wk / kf * P * th1^2;
      Ci[n] := Ci[n] + wk / kf * P * thI^2;
    end for;
  end for;
  for n in 1:nl loop
    D := 4.0 * GeneralRelativity.Constants.pi * ells[n] * (ells[n] + 1.0) / (2.0 * GeneralRelativity.Constants.pi) * (par.TCMB * 1.0e6)^2;
    Dl[n] := D * C0[n];
    Dsw[n] := D * Csw[n];
    Ddop[n] := D * Cd[n];
    Disw[n] := D * Ci[n];
  end for;
  annotation(Documentation(info="<html>
<p>The temperature multipoles of a mode of wavenumber <code>k</code> seen today are</p>
<p><code>&Theta;<sub>l</sub>(k) = &int;d&eta; { g [ (&Theta;<sub>0</sub> + &psi;) j<sub>l</sub>(x) + (&theta;<sub>b</sub>/k) j<sub>l</sub>'(x) ] + e<sup>&minus;&tau;</sup>(&psi;' + &phi;') j<sub>l</sub>(x) }</code>,
<code>x = k(&eta;<sub>0</sub> &minus; &eta;)</code>,</p>
<p>plus the polarization term <code>g&Pi;(j<sub>l</sub> + 3j<sub>l</sub>'')/16</code>
(<code>&Pi; = F<sub>&gamma;2</sub> + G<sub>&gamma;0</sub> + G<sub>&gamma;2</sub></code>): the line-of-sight integral of
Seljak &amp; Zaldarriaga (1996). The terms are the photon's temperature and gravitational redshift where it last
scattered (weighted by the visibility <code>g</code>), the Doppler shift of the moving plasma, the anisotropy of
Thomson scattering, and the integrated Sachs-Wolfe effect of decaying potentials on the way to us. Diffusion
(Silk) damping is in the sources, from the photon Boltzmann hierarchy. Then</p>
<p><code>C<sub>l</sub> = 4&pi; &int; dk/k &Delta;<sub>&Rscr;</sub>&sup2;(k) |&Theta;<sub>l</sub>(k)|&sup2;</code>,
<code>&Delta;<sub>&Rscr;</sub>&sup2; = A<sub>s</sub>(k/k<sub>pivot</sub>)<sup>n<sub>s</sub>&minus;1</sup></code>.</p>
<p>The sources (from <a href=\"modelica://GeneralRelativity.Cosmology.acousticSources\">acousticSources</a>) are
smooth in <code>k</code> and are spline-interpolated to a fine grid that resolves the oscillation of
<code>j<sub>l</sub>(k(&eta;<sub>0</sub> &minus; &eta;<sub>*</sub>))</code> (period about <code>2&pi;/14000</code> per
Mpc). The spherical Bessel functions come from a table on a uniform <code>x</code> grid (cubic Hermite
interpolation; <code>j<sub>l</sub>''</code> from Bessel's equation). The single-term spectra are each source
alone: the total also contains their cross terms, so they do not add up to it. Reionization only enters
through <code>e<sup>&minus;&tau;</sup></code> in the visibility (the suppression <code>e<sup>&minus;2&tau;<sub>re</sub></sup></code>);
the power it regenerates at <code>l &lt; 20</code> is not included.</p>
</html>"));
end cmbSpectrum;
