// The cosmic microwave background (Examples.CosmicMicrowaveBackground): the temperature power spectrum of
// the model against the Planck 2018 data, a sky map drawn from the model spectrum, the sound waves of the
// photon-baryon plasma before recombination, and recombination itself.
//
//   power spectrum: D_l = l(l+1) C_l / 2 pi of the model (and its Sachs-Wolfe, Doppler and ISW parts), the
//     Planck binned TT data with error bars, the Planck best-fit LCDM curve, peaks and residuals;
//   sky: a Gaussian random sky with the model's C_l (seeded), synthesized from spherical harmonics up to
//     l_max, optionally restricted to a range of l to show which scales make which spots; on a Mollweide
//     map or a sphere;
//   sound waves: the radial profile r^2 delta(r) of dark matter, baryons, photons and neutrinos around a
//     point-like initial perturbation (Eisenstein, Seo & White 2007), the Fourier transform of the model's
//     transfer functions at each snapshot: the plasma's sound shell moves out until recombination and stalls
//     near the sound horizon, ~150 Mpc, where dark matter then also collects; and the same in k;
//   recombination: free-electron fraction, optical depth and visibility function against redshift.
(function () {
  // Planck-like diverging colour map for temperature (-1..1).
  const PLANCK = [[-1, [0.0, 0.0, 0.55]], [-0.55, [0.15, 0.45, 1.0]], [-0.15, [0.75, 0.9, 1.0]], [0, [1.0, 0.98, 0.9]],
    [0.25, [1.0, 0.8, 0.45]], [0.6, [1.0, 0.45, 0.1]], [1, [0.55, 0.0, 0.0]]];
  function planckColour(v) {
    v = Math.max(-1, Math.min(1, v));
    for (let i = 1; i < PLANCK.length; i++) {
      if (v <= PLANCK[i][0]) {
        const [x0, c0] = PLANCK[i - 1], [x1, c1] = PLANCK[i], f = (v - x0) / (x1 - x0);
        return [c0[0] + f * (c1[0] - c0[0]), c0[1] + f * (c1[1] - c0[1]), c0[2] + f * (c1[2] - c0[2])];
      }
    }
    return PLANCK[PLANCK.length - 1][1];
  }
  const interp = (xs, ys, x) => {
    if (x <= xs[0]) return ys[0];
    if (x >= xs[xs.length - 1]) return ys[ys.length - 1];
    let lo = 0, hi = xs.length - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (xs[mid] <= x) lo = mid; else hi = mid; }
    return ys[lo] + (x - xs[lo]) / (xs[hi] - xs[lo]) * (ys[hi] - ys[lo]);
  };

  // In-place radix-2 FFT (sign +1: synthesis sum_m X_m e^{+2 pi i m j / N}).
  function fft(re, im) {
    const n = re.length;
    for (let i = 1, j = 0; i < n; i++) {
      let bit = n >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit;
      if (i < j) { [re[i], re[j]] = [re[j], re[i]]; [im[i], im[j]] = [im[j], im[i]]; }
    }
    for (let len = 2; len <= n; len <<= 1) {
      const ang = 2 * Math.PI / len, wr = Math.cos(ang), wi = Math.sin(ang);
      for (let i = 0; i < n; i += len) {
        let cr = 1, ci = 0;
        for (let k = 0; k < len / 2; k++) {
          const a = i + k, b = a + len / 2;
          const tr = re[b] * cr - im[b] * ci, ti = re[b] * ci + im[b] * cr;
          re[b] = re[a] - tr; im[b] = im[a] - ti; re[a] += tr; im[a] += ti;
          const nr = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = nr;
        }
      }
    }
  }

  // A Gaussian random sky with power C[l] (uK^2, l = 0..lmax; zero outside [lmin, lmax]) on an equirectangular
  // grid: nTheta rows (colatitude cell centres) by nPhi columns. Orthonormal spherical harmonics from the
  // standard stable recurrences for the normalized associated Legendre functions.
  function skyMap(C, lmin, lmax, seed) {
    let s = seed >>> 0 || 1;
    const rnd = () => { s = (s + 0x6D2B79F5) >>> 0; let t = s; t = Math.imul(t ^ (t >>> 15), t | 1); t ^= t + Math.imul(t ^ (t >>> 7), t | 61); return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
    const gauss = () => Math.sqrt(-2 * Math.log(rnd() + 1e-300)) * Math.cos(2 * Math.PI * rnd());
    const L = lmax + 1;
    const are = new Float64Array(L * L), aim = new Float64Array(L * L);
    for (let l = 2; l <= lmax; l++) {
      const sd = l >= lmin ? Math.sqrt(Math.max(C[l], 0)) : 0;
      are[l * L] = sd * gauss();
      for (let m = 1; m <= l; m++) { are[l * L + m] = sd * gauss() / Math.SQRT2; aim[l * L + m] = sd * gauss() / Math.SQRT2; }
    }
    // Recurrence coefficients a_lm = sqrt((4 l^2 - 1)/(l^2 - m^2)).
    const coef = new Float64Array(L * L);
    for (let m = 0; m <= lmax; m++) for (let l = m + 1; l <= lmax; l++) coef[l * L + m] = Math.sqrt((4 * l * l - 1) / (l * l - m * m));
    const nTheta = Math.max(64, lmax), nPhi = 1 << Math.ceil(Math.log2(2 * lmax + 2));
    const map = new Float32Array(nTheta * nPhi);
    const re = new Float64Array(nPhi), im = new Float64Array(nPhi);
    for (let i = 0; i < nTheta; i++) {
      const th = Math.PI * (i + 0.5) / nTheta, x = Math.cos(th), sn = Math.sin(th);
      re.fill(0); im.fill(0);
      let pmm = Math.sqrt(1 / (4 * Math.PI));
      for (let m = 0; m <= lmax; m++) {
        if (m > 0) pmm *= -Math.sqrt((2 * m + 1) / (2 * m)) * sn;
        if (Math.abs(pmm) < 1e-280) break;
        let br = are[m * L + m] * pmm, bi = aim[m * L + m] * pmm;
        let p0 = pmm, p1 = m < lmax ? Math.sqrt(2 * m + 3) * x * pmm : 0;
        if (m < lmax) { br += are[(m + 1) * L + m] * p1; bi += aim[(m + 1) * L + m] * p1; }
        for (let l = m + 2; l <= lmax; l++) {
          const p2 = coef[l * L + m] * (x * p1 - p0 / coef[(l - 1) * L + m]);
          br += are[l * L + m] * p2; bi += aim[l * L + m] * p2;
          p0 = p1; p1 = p2;
        }
        const w = m === 0 ? 1 : 2;
        re[m] = w * br; im[m] = w * bi;
      }
      fft(re, im);
      for (let j = 0; j < nPhi; j++) map[i * nPhi + j] = re[j];
    }
    let sum2 = 0;
    for (let k = 0; k < map.length; k++) sum2 += map[k] * map[k];
    return { map, nTheta, nPhi, rms: Math.sqrt(sum2 / map.length) };
  }

  (window.GRStyles = window.GRStyles || []).push({
    id: 'cmb',
    label: 'Cosmic microwave background',
    order: 10,
    description: 'The CMB from the model: temperature power spectrum against the Planck 2018 data, a random sky with that spectrum, the sound waves in the photon-baryon plasma before recombination, and recombination itself.',
    options: [
      { id: 'panel', label: 'view', type: 'select', value: 0, choices: [
        { value: 0, label: 'power spectrum vs Planck' }, { value: 1, label: 'sky map (Mollweide)' }, { value: 2, label: 'sky on a sphere' },
        { value: 3, label: 'sound waves' }, { value: 4, label: 'recombination' }] },
      { id: 'data', label: 'Planck data', type: 'checkbox', value: 1 },
      { id: 'theory', label: 'Planck best fit', type: 'checkbox', value: 1 },
      { id: 'components', label: 'SW / Doppler / ISW', type: 'checkbox', value: 0 },
      { id: 'skyL', label: 'sky l max', type: 'range', min: 64, max: 1024, step: 32, value: 512, help: 'Highest multipole in the sky map (finer spots, slower)' },
      { id: 'lMin', label: 'sky l from', type: 'range', min: 2, max: 1000, step: 1, value: 2, help: 'Keep only multipoles from here ...' },
      { id: 'lMax', label: 'sky l to', type: 'range', min: 2, max: 1024, step: 1, value: 1024, help: '... to here' },
      { id: 'seed', label: 'sky seed', type: 'range', min: 1, max: 100, step: 1, value: 7 },
      { id: 'rMax', label: 'radius (Mpc)', type: 'range', min: 60, max: 400, step: 10, value: 220, help: 'Sound waves: radial range' },
    ],
    available(data) {
      return data.records.cmbSpectrum ? null : 'needs a cmbSpectrum record (Examples.CosmicMicrowaveBackground)';
    },
    create(ctx) {
      const { THREE, data, options } = ctx;
      const rec = (name) => (data.records[name] || [])[0];
      const spec = rec('cmbSpectrum'), dat = rec('cmbData'), theory = rec('cmbTheory'), cmp = rec('cmbComparison') || {};
      const cos = rec('cosmology') || {}, recomb = rec('recombination'), ac = rec('acoustic');
      const layer = GRPlot.layer(ctx);
      const group = new THREE.Group();
      const Tcmb = cos.Tcmb || 2.7255;

      // Model C_l at every integer l (uK^2) for the sky.
      const lTop = Math.round(Math.max(...spec.ell));
      const Cl = new Float64Array(lTop + 1);
      for (let l = 2; l <= lTop; l++) Cl[l] = 2 * Math.PI * interp(spec.ell, spec.Dl, l) / (l * (l + 1));

      // --- sound waves: transfer functions -> real-space profiles --------------------------------------
      let sound = null;
      if (ac) {
        const k = ac.k, nk = k.length, nt = ac.eta.length;
        const dk = k.map((v, i) => 0.5 * ((k[Math.min(i + 1, nk - 1)] - k[Math.max(i - 1, 0)])));
        // Gaussian smoothing (Mpc): suppresses the k cut-off by e^-18, so it does not ring, yet stays far below the
        // ~150 Mpc sound horizon.
        const smooth = Math.max(6 / k[nk - 1], 8);
        const weight = k.map((v, i) => v * v * dk[i] * Math.exp(-0.5 * (v * smooth) ** 2) / (2 * Math.PI * Math.PI));
        // Sign: dark matter overdense at the centre.
        const centre = ac.deltaC[nt - 1].reduce((s, v, i) => s + v * weight[i], 0);
        const sign = centre < 0 ? -1 : 1;
        // Sound horizon r_s(eta) = int c_s d eta, c_s = 1/sqrt(3(1 + R)), R = 3 rho_b / (4 rho_gamma) ~ a; it stops
        // growing when the baryons are released from the photons (the drag epoch, z_drag).
        const omegaGamma = (cos.omegaR || 9.1e-5) / (1 + 0.2271 * (cos.Neff || 3.046));
        const Rb = 3 * (cos.omegaB || 0.0493) / (4 * omegaGamma);
        const logEta = ac.eta.map(Math.log), logA = ac.a.map(Math.log);
        const aOfEta = (e) => Math.exp(interp(logEta, logA, Math.log(e)));
        const rs = [];
        let acc = ac.eta[0] / Math.sqrt(3 * (1 + Rb * ac.a[0])), prev = ac.eta[0];
        rs.push(acc);
        for (let i = 1; i < nt; i++) {
          const steps = 40, h = (ac.eta[i] - prev) / steps;
          const aDrag = 1 / (1 + (cos.zDrag || cos.zStar || 1060));
          for (let s = 0; s < steps; s++) { const e = prev + (s + 0.5) * h, a = aOfEta(e); if (a < aDrag) acc += h / Math.sqrt(3 * (1 + Rb * a)); }
          prev = ac.eta[i];
          rs.push(acc);
        }
        sound = { k, nk, nt, weight, sign, rs, smooth };
      }
      const profileCache = new Map();
      function profiles(rMax) {
        const key = rMax;
        if (profileCache.has(key)) return profileCache.get(key);
        const nr = 240, r = Array.from({ length: nr }, (_, i) => rMax * (i + 0.5) / nr);
        const j0 = r.map((rr) => sound.k.map((kk) => { const x = kk * rr; return x < 1e-6 ? 1 : Math.sin(x) / x; }));
        const out = { r, j0 };
        profileCache.set(key, out);
        return out;
      }
      const transform = (P, T) => P.j0.map((row) => row.reduce((s, v, i) => s + v * T[i] * sound.weight[i], 0) * sound.sign);

      // --- sky ----------------------------------------------------------------------------------------
      let sky = null, skyKey = '', mollImage = null, mollKey = '';
      const sphereMat = new THREE.MeshBasicMaterial({ color: 0xffffff });
      const sphere = new THREE.Mesh(new THREE.SphereGeometry(1, 192, 96), sphereMat);
      group.add(sphere);
      function ensureSky() {
        const lmax = Math.min(Math.round(options.skyL), lTop), lmin = Math.max(2, Math.round(options.lMin)), lmaxF = Math.min(lmax, Math.round(options.lMax));
        const key = [lmax, lmin, lmaxF, options.seed].join('|');
        if (key === skyKey) return;
        skyKey = key;
        const C = new Float64Array(lmax + 1);
        for (let l = lmin; l <= lmaxF; l++) C[l] = Cl[l];
        sky = skyMap(C, lmin, lmax, Math.round(options.seed) * 7919);
        sky.lmin = lmin; sky.lmax = lmaxF;
        // Sphere texture (row 0 = north pole at the top).
        const { nTheta, nPhi, map, rms } = sky;
        const rgba = new Uint8Array(nTheta * nPhi * 4);
        for (let i = 0; i < nTheta; i++) for (let j = 0; j < nPhi; j++) {
          const c = planckColour(map[i * nPhi + j] / (3 * rms || 1)), o = 4 * ((nTheta - 1 - i) * nPhi + j);
          rgba[o] = 255 * c[0]; rgba[o + 1] = 255 * c[1]; rgba[o + 2] = 255 * c[2]; rgba[o + 3] = 255;
        }
        if (sphereMat.map) sphereMat.map.dispose();
        const tex = new THREE.DataTexture(rgba, nPhi, nTheta, THREE.RGBAFormat);
        tex.colorSpace = THREE.SRGBColorSpace; tex.magFilter = THREE.LinearFilter; tex.minFilter = THREE.LinearFilter; tex.wrapS = THREE.RepeatWrapping;
        tex.needsUpdate = true;
        sphereMat.map = tex; sphereMat.needsUpdate = true;
        mollKey = '';
      }
      // Mollweide image of the sky in the box (cached).
      function mollweide(g, box) {
        const key = skyKey + '|' + [box.x, box.y, box.w, box.h].join(',');
        if (key !== mollKey) {
          mollKey = key;
          const W = Math.round(box.w), H = Math.round(box.h), img = g.createImageData(W, H);
          const { map, nTheta, nPhi, rms } = sky;
          for (let py = 0; py < H; py++) {
            const y = Math.SQRT2 * (1 - 2 * (py + 0.5) / H);
            const aux = Math.asin(Math.max(-1, Math.min(1, y / Math.SQRT2)));
            const lat = Math.asin(Math.max(-1, Math.min(1, (2 * aux + Math.sin(2 * aux)) / Math.PI)));
            const ti = ((Math.PI / 2 - lat) / Math.PI) * nTheta - 0.5;
            for (let px = 0; px < W; px++) {
              const x = 2 * Math.SQRT2 * (2 * (px + 0.5) / W - 1);
              const lon = Math.PI * x / (2 * Math.SQRT2 * Math.cos(aux));
              const o = 4 * (py * W + px);
              if (Math.abs(lon) > Math.PI) { img.data[o + 3] = 0; continue; }
              // Galactic-style orientation: longitude increasing to the left.
              const tj = ((Math.PI - lon) / (2 * Math.PI)) * nPhi;
              const i0 = Math.max(0, Math.min(nTheta - 1, Math.floor(ti))), i1 = Math.min(nTheta - 1, i0 + 1), fi = Math.max(0, Math.min(1, ti - i0));
              const j0 = ((Math.floor(tj) % nPhi) + nPhi) % nPhi, j1 = (j0 + 1) % nPhi, fj = tj - Math.floor(tj);
              const v = (1 - fi) * ((1 - fj) * map[i0 * nPhi + j0] + fj * map[i0 * nPhi + j1]) + fi * ((1 - fj) * map[i1 * nPhi + j0] + fj * map[i1 * nPhi + j1]);
              const c = planckColour(v / (3 * rms || 1));
              img.data[o] = 255 * c[0]; img.data[o + 1] = 255 * c[1]; img.data[o + 2] = 255 * c[2]; img.data[o + 3] = 255;
            }
          }
          mollImage = img;
        }
        g.putImageData(mollImage, Math.round(box.x), Math.round(box.y));
      }
      function colourBar(g, x, y, w, h, lo, hi, label, u) {
        for (let i = 0; i < w; i++) { const c = planckColour(2 * i / (w - 1) - 1); g.fillStyle = `rgb(${255 * c[0]},${255 * c[1]},${255 * c[2]})`; g.fillRect(x + i, y, 1, h); }
        g.fillStyle = '#c9d1e0'; g.font = `${Math.round(11 * u)}px system-ui`; g.textBaseline = 'top';
        g.textAlign = 'left'; g.fillText(`${Math.round(lo)} μK`, x, y + h + 3 * u);
        g.textAlign = 'right'; g.fillText(`+${Math.round(hi)} μK`, x + w, y + h + 3 * u);
        g.textAlign = 'center'; g.fillText(label, x + w / 2, y + h + 3 * u);
      }

      // --- l axis like Planck's: logarithmic up to 30, linear above -------------------------------------
      const lHi = 2500;
      const lMap = (l) => (l < 30 ? 0.14 * Math.log(Math.max(l, 2) / 2) / Math.log(15) : 0.14 + 0.86 * (l - 30) / (lHi - 30));
      const lTicks = [2, 10, 30, 500, 1000, 1500, 2000, 2500];

      let lastKey = '';
      function draw(r) {
        const { w, h, u } = layer.begin(r);
        const panel = Math.round(options.panel);
        const key = [w, h, panel, panel === 3 ? ctx.clock : '', JSON.stringify(options)].join('|');
        if (key === lastKey && panel !== 2) return;
        lastKey = key;
        const g = layer.g;
        if (panel === 2) g.clearRect(0, 0, w, h); else { g.fillStyle = '#0b0e14'; g.fillRect(0, 0, w, h); }
        g.save();
        g.font = `${Math.round(13 * u)}px system-ui`; g.fillStyle = '#e6e9ef'; g.textBaseline = 'top';
        g.fillText(['Temperature power spectrum: model vs Planck 2018', 'A sky with the model spectrum (random realization)', '', 'Sound waves in the early universe', 'Recombination'][panel], 70 * u, 8 * u);
        g.restore();
        if (panel === 0) drawSpectrum(g, w, h, u);
        else if (panel === 1 || panel === 2) drawSky(g, w, h, u, panel);
        else if (panel === 3) drawSound(g, w, h, u);
        else drawRecombination(g, w, h, u);
      }

      function drawSpectrum(g, w, h, u) {
        const box = { x: 70 * u, y: 34 * u, w: w - 100 * u, h: 0.66 * h - 40 * u };
        const box2 = { x: box.x, y: box.y + box.h + 36 * u, w: box.w, h: h - box.y - box.h - 80 * u };
        const yTop = 6500;
        const m = GRPlot.axes(g, box, { lo: 2, hi: lHi, map: lMap, ticks: lTicks, noLabels: true }, { lo: 0, hi: yTop, label: 'D_l = l(l+1) C_l / 2π  (μK²)' }, u);
        if (options.theory && theory) GRPlot.line(g, m, theory.ell, theory.Dl, { colour: '#8a93a6', width: 1.2 * u, dash: [6 * u, 4 * u] }, box);
        if (options.components) {
          [['sw', '#7ee0a1'], ['doppler', '#c792ea'], ['isw', '#ff6b8a']].forEach(([key, colour]) => {
            if (spec[key]) GRPlot.line(g, m, spec.ell, spec[key], { colour, width: 1.2 * u }, box);
          });
        }
        GRPlot.line(g, m, spec.ell, spec.Dl, { colour: '#ffb454', width: 2.4 * u }, box);
        if (options.data && dat) {
          g.save();
          g.beginPath(); g.rect(box.x, box.y, box.w, box.h); g.clip();
          for (let i = 0; i < dat.ell.length; i++) {
            const x = m.X(dat.ell[i]), y = m.Y(dat.Dl[i]);
            g.strokeStyle = 'rgba(230,233,239,0.7)'; g.lineWidth = Math.max(1, u);
            g.beginPath(); g.moveTo(x, m.Y(dat.Dl[i] - dat.errLo[i])); g.lineTo(x, m.Y(dat.Dl[i] + dat.errHi[i])); g.stroke();
            g.fillStyle = '#e6e9ef'; g.beginPath(); g.arc(x, y, 2 * u, 0, 2 * Math.PI); g.fill();
          }
          g.restore();
        }
        // Peaks.
        g.save();
        g.font = `${Math.round(10 * u)}px system-ui`; g.textAlign = 'center';
        (cmp.peaksModel || []).forEach(([l, D], i) => {
          g.fillStyle = '#ffb454'; g.fillText(`${Math.round(l)}`, m.X(l), m.Y(D) - 16 * u);
          const d = (cmp.peaksData || [])[i];
          if (d) { g.fillStyle = '#e6e9ef'; g.fillText(`(${Math.round(d[0])})`, m.X(l), m.Y(D) - 28 * u); }
        });
        g.restore();
        const legend = [{ label: 'this model (Modelica)', colour: '#ffb454' }];
        if (options.theory && theory) legend.push({ label: 'Planck 2018 best-fit LCDM', colour: '#8a93a6', dash: [6 * u, 4 * u] });
        if (options.data && dat) legend.push({ label: 'Planck 2018 data (binned above l = 30)', colour: '#e6e9ef', dot: true });
        if (options.components) legend.push({ label: 'Sachs-Wolfe + monopole', colour: '#7ee0a1' }, { label: 'Doppler', colour: '#c792ea' }, { label: 'integrated Sachs-Wolfe', colour: '#ff6b8a' });
        GRPlot.legend(g, box.x + box.w - 300 * u, box.y + 16 * u, legend, u);
        // Residuals (data - model) / error.
        if (dat) {
          const mr = GRPlot.axes(g, box2, { lo: 2, hi: lHi, map: lMap, ticks: lTicks, label: 'multipole l' }, { lo: -5, hi: 5, ticks: [-4, -2, 0, 2, 4], label: '(data − model)/σ' }, u);
          g.save();
          g.fillStyle = 'rgba(95,179,255,0.10)'; g.fillRect(box2.x, mr.Y(2), box2.w, mr.Y(-2) - mr.Y(2));
          g.fillStyle = 'rgba(95,179,255,0.14)'; g.fillRect(box2.x, mr.Y(1), box2.w, mr.Y(-1) - mr.Y(1));
          g.beginPath(); g.rect(box2.x, box2.y, box2.w, box2.h); g.clip();
          for (let i = 0; i < dat.ell.length; i++) {
            const model = interp(spec.ell, spec.Dl, dat.ell[i]), d = dat.Dl[i] - model, sig = d > 0 ? dat.errLo[i] : dat.errHi[i];
            const z = Math.max(-5, Math.min(5, d / Math.max(sig, 1e-9)));
            g.fillStyle = Math.abs(z) < 2 ? '#e6e9ef' : '#ff6b8a';
            g.beginPath(); g.arc(mr.X(dat.ell[i]), mr.Y(z), 2 * u, 0, 2 * Math.PI); g.fill();
          }
          g.restore();
        }
        // Summary.
        g.save();
        g.font = `${Math.round(11 * u)}px system-ui`; g.fillStyle = '#c9d1e0'; g.textBaseline = 'top';
        const lines = [];
        if (cmp.chi2 !== undefined) lines.push(`χ² = ${cmp.chi2.toFixed(0)} for ${cmp.nBins} Planck bins (diagonal errors)`);
        if (cmp.rmsVsTheory !== undefined) lines.push(`rms difference from the Planck best fit: ${(100 * cmp.rmsVsTheory).toFixed(1)}%`);
        if (cos.thetaStar) lines.push(`100 θ* = ${cos.thetaStar.toFixed(4)} (Planck 1.0411), z* = ${ctx.fmt(cos.zStar)}, r_s = ${ctx.fmt(cos.rsStar)} Mpc`);
        lines.forEach((t, i) => g.fillText(t, box.x + box.w - 300 * u, box.y + 24 * u + 16 * u * legend.length + 15 * u * i));
        g.restore();
      }

      function drawSky(g, w, h, u, panel) {
        ensureSky();
        const W = Math.min(w - 40 * u, 2 * (h - 110 * u)), H = W / 2;
        const box = { x: (w - W) / 2, y: 40 * u, w: W, h: H };
        if (panel === 1) mollweide(g, box);
        const sig = 3 * sky.rms;
        colourBar(g, w / 2 - 160 * u, (panel === 1 ? box.y + H : h - 70 * u) + 14 * u, 320 * u, 10 * u, -sig, sig, 'ΔT', u);
        g.save();
        g.font = `${Math.round(11 * u)}px system-ui`; g.fillStyle = '#c9d1e0'; g.textAlign = 'center'; g.textBaseline = 'top';
        g.fillText(`Multipoles ${sky.lmin} to ${sky.lmax} of the model spectrum, seed ${Math.round(options.seed)}; rms ${sky.rms.toFixed(0)} μK. `
          + (sky.lmin < 300 && sky.lmax > 150 ? 'The typical ~1° spots are the first acoustic peak (l ≈ 220).' : 'Change the l range to see which scales make which spots.'), w / 2, (panel === 1 ? box.y + H : h - 70 * u) + 44 * u);
        g.restore();
      }

      function drawSound(g, w, h, u) {
        if (!sound) { g.fillStyle = '#c9d1e0'; g.fillText('The scene has no acoustic record.', 70 * u, 60 * u); return; }
        const f = Math.max(0, Math.min(sound.nt - 1, ctx.clock)), i0 = Math.min(sound.nt - 2, Math.floor(f)), t = f - i0;
        const mix = (arr) => arr[i0].map((v, k) => (1 - t) * v + t * arr[i0 + 1][k]);
        const species = [['dark matter', '#5fb3ff', ac.deltaC], ['baryons', '#ffb454', ac.deltaB], ['photons', '#ff6b8a', ac.deltaG], ['neutrinos', '#7ee0a1', ac.deltaNu]].filter((s) => s[2]);
        const z = Math.exp((1 - t) * Math.log(ac.z[i0] + 1) + t * Math.log(ac.z[i0 + 1] + 1)) - 1;
        const rs = (1 - t) * sound.rs[i0] + t * sound.rs[i0 + 1];
        const P = profiles(options.rMax);
        const curves = species.map(([label, colour, arr]) => ({ label, colour, T: mix(arr) }));
        const rp = curves.map((c) => transform(P, c.T).map((d, i) => d * P.r[i] * P.r[i]));
        const top = Math.max(1e-30, ...rp.flat().map(Math.abs)) * 1.1;
        const box = { x: 70 * u, y: 34 * u, w: w - 100 * u, h: 0.55 * h - 40 * u };
        const m = GRPlot.axes(g, box, { lo: 0, hi: options.rMax, label: 'comoving radius r (Mpc)' }, { lo: -top, hi: top, label: 'r² δ(r)  (scaled to fit)', format: () => '' }, u);
        curves.forEach((c, k) => GRPlot.line(g, m, P.r, rp[k], { colour: c.colour, width: 2 * u }, box));
        g.save();
        g.strokeStyle = 'rgba(255,255,255,0.5)'; g.setLineDash([5 * u, 4 * u]); g.lineWidth = 1.2 * u;
        g.beginPath(); g.moveTo(m.X(rs), box.y); g.lineTo(m.X(rs), box.y + box.h); g.stroke();
        g.fillStyle = '#c9d1e0'; g.font = `${Math.round(11 * u)}px system-ui`;
        g.textAlign = 'right';
        g.fillText(`sound horizon ${rs.toFixed(0)} Mpc${z < (cos.zDrag || 1060) ? ' (frozen at the drag epoch)' : ''}`, m.X(rs) - 4 * u, box.y + box.h - 14 * u);
        g.restore();
        GRPlot.legend(g, box.x + box.w - 160 * u, box.y + 16 * u, curves.map((c) => ({ label: c.label, colour: c.colour })), u);
        // k space below.
        const box2 = { x: box.x, y: box.y + box.h + 50 * u, w: box.w, h: h - box.y - box.h - 100 * u };
        // Scaled to the oscillating species; dark matter, which keeps growing, may run off the top.
        const osc = curves.filter((c) => c.label !== 'dark matter');
        const kTop = Math.max(1e-30, ...(osc.length ? osc : curves).flatMap((c) => c.T.map(Math.abs))) * 1.3;
        const mk = GRPlot.axes(g, box2, { lo: 0, hi: sound.k[sound.nk - 1], label: 'wavenumber k (1/Mpc)' }, { lo: -kTop, hi: kTop, label: 'δ(k), per unit curvature' }, u);
        curves.forEach((c) => GRPlot.line(g, mk, sound.k, c.T.map((v) => v * sound.sign), { colour: c.colour, width: 1.6 * u }, box2));
        g.save();
        g.font = `${Math.round(12 * u)}px system-ui`; g.fillStyle = '#e6e9ef'; g.textBaseline = 'top';
        g.fillText(`z = ${ctx.fmt(z)}  (${z > (cos.zStar || 1090) ? 'plasma: photons and baryons oscillate together' : 'after recombination: baryons fall into the dark matter wells'})`, box.x + 10 * u, box.y + 10 * u);
        g.font = `${Math.round(10 * u)}px system-ui`; g.fillStyle = '#8a93a6';
        g.fillText(`Profiles around a point-like initial perturbation (Fourier transform of the transfer functions, smoothed over ${sound.smooth.toFixed(0)} Mpc).`, box.x + 10 * u, box.y + 28 * u);
        g.restore();
      }

      function drawRecombination(g, w, h, u) {
        if (!recomb) { g.fillStyle = '#c9d1e0'; g.fillText('The scene has no recombination record.', 70 * u, 60 * u); return; }
        const box = { x: 70 * u, y: 34 * u, w: w - 140 * u, h: h - 100 * u };
        const zs = recomb.z.map((z) => z + 1);
        const m = GRPlot.axes(g, box, { lo: 1, hi: Math.max(...zs), log: true, label: '1 + z' }, { lo: 1e-4, hi: 2, log: true, label: 'free electrons per hydrogen nucleus  x_e' }, u);
        GRPlot.line(g, m, zs, recomb.xe.map((x) => Math.max(x, 1e-4)), { colour: '#ffb454', width: 2.2 * u }, box);
        if (recomb.tau) GRPlot.line(g, m, zs, recomb.tau.map((x) => Math.max(x, 1e-4)), { colour: '#5fb3ff', width: 1.8 * u }, box);
        if (recomb.g) {
          const gMax = Math.max(...recomb.g);
          const mg = { X: m.X, Y: (v) => box.y + box.h * (1 - v / (gMax * 1.1)) };
          GRPlot.line(g, mg, zs, recomb.g, { colour: '#7ee0a1', width: 2 * u }, box);
        }
        if (cos.zStar) {
          g.save(); g.strokeStyle = 'rgba(255,255,255,0.4)'; g.setLineDash([4 * u, 4 * u]);
          g.beginPath(); g.moveTo(m.X(cos.zStar + 1), box.y); g.lineTo(m.X(cos.zStar + 1), box.y + box.h); g.stroke();
          g.fillStyle = '#c9d1e0'; g.font = `${Math.round(11 * u)}px system-ui`; g.fillText(`z* = ${ctx.fmt(cos.zStar)}`, m.X(cos.zStar + 1) + 4 * u, box.y + 14 * u);
          g.restore();
        }
        GRPlot.legend(g, box.x + 12 * u, box.y + 16 * u, [{ label: 'free-electron fraction x_e', colour: '#ffb454' }, { label: 'optical depth τ to today', colour: '#5fb3ff' }, { label: 'visibility g (where CMB photons last scattered; linear scale)', colour: '#7ee0a1' }], u);
      }

      if (ac) ctx.setTimeRange(0, ac.eta.length - 1);
      const notes = [];
      if (cmp.chi2 !== undefined) notes.push(`Model vs Planck 2018 binned TT: χ² = ${cmp.chi2.toFixed(1)} for ${cmp.nBins} bins (diagonal errors only, so indicative).`);
      notes.push('Data: Planck Collaboration 2020, A&A 641, A1/A6 (Planck Legacy Archive, PR3).');
      if (ac) notes.push('The clock steps through the sound-wave snapshots (sound waves view).');
      ctx.setNotes(notes.join('\n'));
      return {
        group,
        hidesBodies: true,
        hidesTracks: true,
        ownClock: !!ac,
        clockLabel: 'step',
        place(x) { return ctx.toWorld(x[1], x[2], x[3]); },
        camera() { return { position: new THREE.Vector3(0, 0.6, 2.6), target: new THREE.Vector3(0, 0, 0) }; },
        render(r, camera, scene3) {
          const panel = Math.round(options.panel);
          sphere.visible = panel === 2;
          if (panel === 2) {
            ensureSky();
            r.render(scene3, camera);
            draw(r);
            layer.render(r, true);
          } else {
            draw(r);
            layer.render(r);
          }
        },
        readout() {
          const lines = [];
          if (cmp.peaksModel) lines.push('peaks (l, D_l): model ' + cmp.peaksModel.map(([l, D]) => `${Math.round(l)}/${Math.round(D)}`).join(', '));
          if (cmp.peaksData) lines.push('               Planck ' + cmp.peaksData.map(([l, D]) => `${Math.round(l)}/${Math.round(D)}`).join(', '));
          if (sound && Math.round(options.panel) === 3) {
            const i = Math.max(0, Math.min(sound.nt - 1, Math.round(ctx.clock)));
            lines.push(`snapshot ${i + 1} of ${sound.nt}: z = ${ctx.fmt(ac.z[i])}, eta = ${ctx.fmt(ac.eta[i])} Mpc`);
          }
          return lines;
        },
        onOption() { lastKey = ''; return true; },
        dispose() { layer.dispose(); if (sphereMat.map) sphereMat.map.dispose(); },
      };
    },
  });
})();
