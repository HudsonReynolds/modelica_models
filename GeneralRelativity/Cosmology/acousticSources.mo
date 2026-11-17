within GeneralRelativity.Cosmology;
function acousticSources "Evolve the perturbations of every wavenumber: CMB sources at recombination, ISW source, transfer-function snapshots"
  input Parameters par;
  input Real k[:](each unit = "1/Mpc") "Wavenumbers";
  input Real etaRec[:] "Increasing conformal times spanning last scattering";
  input Real etaIsw[:] "Increasing conformal times after last scattering, to today";
  input Real etaSnap[:] "Increasing conformal times of the transfer-function snapshots";
  input Real etaOff "Earliest time at which radiation perturbations may be dropped (a few times eta*)";
  input Real eta0 "Conformal time today";
  input Real lnEtaStart "Thermal tables: ln(eta) of the first entry";
  input Real dLnEta "Thermal tables: spacing in ln(eta)";
  input Real dtauTable[:] "Thomson scattering rate, 1/Mpc";
  input Real kDTable[size(dtauTable, 1)] "Diffusion-damping wavenumber, 1/Mpc";
  input Real aTable[size(dtauTable, 1)] "Scale factor";
  input Real etaTightMax = 1.0e30 "Tight coupling ends here at the latest";
  output Real rec[size(k, 1), size(etaRec, 1), 5] "delta_gamma/4, psi, theta_b/k, psi' + phi', Pi";
  output Real isw[size(k, 1), size(etaIsw, 1)] "psi' + phi'";
  output Real snap[size(k, 1), size(etaSnap, 1), 5] "delta_c, delta_b, delta_gamma, delta_nu, phi";
protected
  Real H0, Og, On, Ob, Oc, OL;
  Real r[size(etaRec, 1), 5];
  Real s[size(etaSnap, 1), 5];
  Real w[size(etaIsw, 1)];
  Real etaTight;
algorithm
  (H0, Og, On, Ob, Oc, OL) := densities(par);
  for i in 1:size(k, 1) loop
    (r, w, s, etaTight) := evolveMode(k[i], etaRec, etaIsw, etaSnap, etaOff, eta0, lnEtaStart, dLnEta, dtauTable, kDTable, aTable,
      H0, Og, On, Ob, Oc, par.omegaK, OL, etaTightMax = etaTightMax);
    rec[i, :, :] := r;
    isw[i, :] := w;
    snap[i, :, :] := s;
  end for;
end acousticSources;
