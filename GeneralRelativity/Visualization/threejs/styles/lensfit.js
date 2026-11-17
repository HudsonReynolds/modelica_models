// Lens-fit style: the result of a lens inversion (GeneralRelativity.Lensing.Inversion). Panels show the observed
// image, the simulation of the first radius-based estimate, the fitted simulation and the residual in units of the
// pixel noise; below them the cost and the mass estimate per iteration, against the truth when the scene knows it.
// Records: lensImage (name, u, v, values[i][j] at (u[i], v[j]) in arcseconds) and lensFit (parameters, errors, history).
(function () {
  const VERTEX = 'varying vec2 vUv; void main() { vUv = uv; gl_Position = vec4(position.xy, 0.0, 1.0); }';
  const FRAGMENT = /* glsl */`
    varying vec2 vUv;
    uniform sampler2D uTex;
    uniform int uMode;            // 0 brightness, 1 residual
    uniform float uExposure, uPeak, uRange, uLog;
    uniform vec3 uRing;           // centre x, y and radius of the fitted Einstein ring, in texture units (0 radius: none)
    uniform vec2 uSpan;           // extent of the panel in arcsec per texture unit, for the ring's aspect
    vec3 heat(float t) {
      t = clamp(t, 0.0, 1.0);
      vec3 c = t < 0.33 ? mix(vec3(0.02, 0.01, 0.08), vec3(0.45, 0.05, 0.45), t / 0.33)
             : t < 0.66 ? mix(vec3(0.45, 0.05, 0.45), vec3(0.98, 0.45, 0.1), (t - 0.33) / 0.33)
                        : mix(vec3(0.98, 0.45, 0.1), vec3(1.0, 1.0, 0.85), (t - 0.66) / 0.34);
      return c;
    }
    void main() {
      float x = texture2D(uTex, vUv).r;
      vec3 col;
      if (uMode == 0) {
        float s = x / max(uPeak, 1e-30) * uExposure;
        float t = uLog > 0.5 ? log(1.0 + 30.0 * max(s, 0.0)) / log(31.0) : s;
        col = heat(t);
      } else {
        float t = clamp(x / uRange, -1.0, 1.0);
        vec3 mid = vec3(0.13, 0.14, 0.17);
        col = t < 0.0 ? mix(mid, vec3(0.25, 0.55, 1.0), -t) : mix(mid, vec3(1.0, 0.35, 0.2), t);
      }
      if (uRing.z > 0.0) {
        vec2 d = (vUv - uRing.xy) * uSpan;
        float r = length(d);
        float w = fwidth(r) * 1.2;
        col = mix(col, vec3(0.4, 0.95, 1.0), (1.0 - smoothstep(0.0, w, abs(r - uRing.z))) * 0.8);
      }
      gl_FragColor = vec4(col, 1.0);
    }`;
  const ORDER = ['observed', 'initial', 'model', 'residual'];
  const TITLES = { observed: 'observed image', initial: 'first simulation (radius estimate)', model: 'fitted simulation', residual: 'residual (observed - fitted)' };

  window.GRStyles = window.GRStyles || [];
  window.GRStyles.push({
    id: 'lensfit',
    label: 'Lens inversion',
    order: 8,
    description: 'The mass behind a lensed image: the observed picture, the simulation of the first estimate from the ring radius, the simulation after the least-squares fit and what is left over, with the cost and mass per iteration. Scroll to change the exposure.',
    options: [
      { id: 'panels', label: 'panels', type: 'select', choices: [{ value: 0, label: 'all' }, { value: 1, label: 'observed' }, { value: 2, label: 'first simulation' }, { value: 3, label: 'fitted' }, { value: 4, label: 'residual' }, { value: 5, label: 'observed | fitted' }], value: 0 },
      { id: 'exposure', label: 'exposure', type: 'range', min: 0.2, max: 8, step: 0.1, value: 1.5 },
      { id: 'log', label: 'log stretch', type: 'checkbox', value: 1 },
      { id: 'range', label: 'residual range (sigma)', type: 'range', min: 1, max: 20, step: 0.5, value: 5 },
      { id: 'ring', label: 'fitted Einstein ring', type: 'checkbox', value: 1 },
      { id: 'history', label: 'iteration plots', type: 'checkbox', value: 1 },
      { id: 'labels', label: 'labels', type: 'checkbox', value: 1 },
    ],
    available(data) {
      return data.records && data.records.lensImage ? null : 'needs lens-inversion images (Examples.LensInversion or Visualization.writeLensImage)';
    },
    create(ctx) {
      const { THREE, data, options, renderer, controls } = ctx;
      controls.enabled = false;
      const previousClear = renderer.getClearColor(new THREE.Color()), previousAlpha = renderer.getClearAlpha();
      const images = {};
      for (const rec of data.records.lensImage) images[rec.name] = rec;
      const fit = (data.records.lensFit || [])[0] || null;
      const names = ORDER.filter(n => images[n]).concat(Object.keys(images).filter(n => !ORDER.includes(n)));
      const quadCamera = new THREE.OrthographicCamera(-1, 1, 1, -1, 0, 1);
      const panels = names.map(name => {
        const rec = images[name];
        const nu = rec.u.length, nv = rec.v.length;
        const tex = new Float32Array(nu * nv * 4);
        let peak = 0;
        for (let i = 0; i < nu; i++) for (let j = 0; j < nv; j++) {
          const val = rec.values[i][j];
          tex[(j * nu + i) * 4] = val;
          if (name !== 'residual' && val > peak) peak = val;
        }
        const texture = new THREE.DataTexture(tex, nu, nv, THREE.RGBAFormat, THREE.FloatType);
        texture.minFilter = texture.magFilter = THREE.NearestFilter;
        texture.needsUpdate = true;
        const du = nu > 1 ? rec.u[1] - rec.u[0] : 1, dv = nv > 1 ? rec.v[1] - rec.v[0] : 1;
        const u0 = rec.u[0] - 0.5 * du, v0 = rec.v[0] - 0.5 * dv, spanU = nu * du, spanV = nv * dv;
        const uniforms = {
          uTex: { value: texture }, uMode: { value: name === 'residual' ? 1 : 0 }, uExposure: { value: 1 }, uPeak: { value: peak || 1 },
          uRange: { value: 5 }, uLog: { value: 1 }, uRing: { value: new THREE.Vector3(0, 0, 0) }, uSpan: { value: new THREE.Vector2(spanU, spanV) },
        };
        const material = new THREE.ShaderMaterial({ vertexShader: VERTEX, fragmentShader: FRAGMENT, uniforms, depthTest: false, depthWrite: false });
        const scene = new THREE.Scene();
        scene.add(new THREE.Mesh(new THREE.PlaneGeometry(2, 2), material));
        return { name, rec, scene, material, texture, uniforms, u0, v0, spanU, spanV, aspect: spanU / spanV, peak, unit: rec.unit || '' };
      });
      // The shared peak keeps the observed and simulated panels on one brightness scale.
      const sharedPeak = Math.max(...panels.filter(p => p.name !== 'residual').map(p => p.peak), 1e-30);
      for (const p of panels) if (p.name !== 'residual') p.uniforms.uPeak.value = sharedPeak;

      // Iteration plots: cost (log scale) and mass estimate against the iteration, with the truth and the first estimate.
      const plotScene = new THREE.Scene();
      const plotCamera = new THREE.OrthographicCamera(0, 1, 1, 0, -1, 1);
      const lineMat = (color, dashed) => dashed ? new THREE.LineDashedMaterial({ color, dashSize: 0.02, gapSize: 0.012 }) : new THREE.LineBasicMaterial({ color });
      const plots = [];
      function addLine(points, color, dashed) {
        const line = new THREE.Line(new THREE.BufferGeometry().setFromPoints(points.map(p => new THREE.Vector3(p[0], p[1], 0))), lineMat(color, dashed));
        if (dashed) line.computeLineDistances();
        plotScene.add(line);
        plots.push(line);
        return line;
      }
      let plotInfo = [];
      if (fit && fit.chi2History && fit.chi2History.length > 1) {
        const hist = fit.chi2History, n = hist.length, last = Math.max(1, fit.iterations || n - 1);
        const used = hist.slice(0, Math.min(n, last + 1));
        const lo = Math.log10(Math.max(Math.min(...used), 1e-30)), hi = Math.log10(Math.max(...used));
        const span = Math.max(hi - lo, 1e-6);
        // Left half: chi2. Right half: mass. Each in a 0..0.46 wide box with a 0.03 margin.
        const box = (k) => ({ x0: 0.03 + 0.5 * k, x1: 0.47 + 0.5 * k, y0: 0.1, y1: 0.92 });
        const b1 = box(0), b2 = box(1);
        for (const b of [b1, b2]) addLine([[b.x0, b.y0], [b.x1, b.y0], [b.x1, b.y1], [b.x0, b.y1], [b.x0, b.y0]], 0x5a6478);
        const sx = (b, i) => b.x0 + (b.x1 - b.x0) * i / Math.max(used.length - 1, 1);
        addLine(used.map((c, i) => [sx(b1, i), b1.y0 + (b1.y1 - b1.y0) * (Math.log10(Math.max(c, 1e-30)) - lo) / span]), 0xffb454);
        if (fit.pixels > 0) {
          const y = b1.y0 + (b1.y1 - b1.y0) * (Math.log10(fit.pixels) - lo) / span;
          if (y > b1.y0 && y < b1.y1) addLine([[b1.x0, y], [b1.x1, y]], 0x8a93a6, true);
        }
        plotInfo.push(`chi2 per iteration (log): ${used[0].toExponential(2)} -> ${used[used.length - 1].toFixed(0)}${fit.pixels ? `, dashed = ${fit.pixels} pixels` : ''}`);
        const mass = (fit.massHistorySolar || []).slice(0, used.length);
        if (mass.length > 1) {
          const all = mass.concat(fit.trueMassSolar > 0 ? [fit.trueMassSolar] : []);
          let mlo = Math.min(...all), mhi = Math.max(...all);
          const pad = 0.1 * Math.max(mhi - mlo, 1e-6 * mhi); mlo -= pad; mhi += pad;
          const my = (m) => b2.y0 + (b2.y1 - b2.y0) * (m - mlo) / (mhi - mlo);
          addLine(mass.map((m, i) => [sx(b2, i), my(m)]), 0x6ad0ff);
          if (fit.trueMassSolar > 0) addLine([[b2.x0, my(fit.trueMassSolar)], [b2.x1, my(fit.trueMassSolar)]], 0x7cff9b, true);
          if (fit.massSigmaSolar > 0) {
            const yl = my(fit.massSolar - fit.massSigmaSolar), yh = my(fit.massSolar + fit.massSigmaSolar);
            addLine([[b2.x1 - 0.01, yl], [b2.x1 - 0.01, yh]], 0xffffff);
          }
          plotInfo.push(`mass inside the Einstein radius per iteration: ${mass[0].toExponential(3)} -> ${mass[mass.length - 1].toExponential(3)} M_sun${fit.trueMassSolar > 0 ? ` (dashed: truth ${fit.trueMassSolar.toExponential(3)})` : ''}`);
        }
      }

      // Labels over the panels, positioned over the canvas.
      const el = renderer.domElement;
      const overlay = document.createElement('div');
      overlay.style.cssText = 'position:fixed;left:0;top:0;pointer-events:none;font:12px/1.3 system-ui,sans-serif;color:#dfe5f0;text-shadow:0 0 4px #000,0 0 2px #000;z-index:5;';
      document.body.appendChild(overlay);
      const labelEls = [];
      function label(k) {
        while (labelEls.length <= k) { const d = document.createElement('div'); d.style.position = 'absolute'; d.style.whiteSpace = 'pre'; overlay.appendChild(d); labelEls.push(d); }
        return labelEls[k];
      }
      const onWheel = (e) => {
        e.preventDefault();
        options.exposure = Math.max(0.2, Math.min(8, options.exposure * Math.exp(-0.0015 * e.deltaY)));
        ctx.refreshReadout();
        clearTimeout(wheelTimer);
        wheelTimer = setTimeout(() => ctx.refreshOptionsUI(), 250);
      };
      let wheelTimer = null;
      el.addEventListener('wheel', onWheel, { passive: false });

      const notes = [];
      if (fit) {
        notes.push(`${fit.model === 'point' ? 'Point-mass' : 'Isothermal-ellipse'} lens fitted in ${fit.iterations} iterations${fit.converged ? '' : ' (iteration limit reached)'}; chi2 ${fit.chi2History[Math.min(fit.iterations, fit.chi2History.length - 1)].toFixed(0)} for ${fit.pixels} pixels, ${(100 * fit.explained).toFixed(1)}% of the image explained, correlation ${fit.correlation.toFixed(4)}.`);
        notes.push(`Mass inside the Einstein radius: ${fit.massInitialSolar.toExponential(3)} M_sun from the ring radius, ${fit.massSolar.toExponential(3)} +/- ${fit.massSigmaSolar.toExponential(2)} after the fit${fit.trueMassSolar > 0 ? `; truth ${fit.trueMassSolar.toExponential(3)} (${((fit.massSolar - fit.trueMassSolar) / fit.massSigmaSolar).toFixed(1)} sigma off)` : ''}.`);
        notes.push('Standard errors come from the curvature of chi2 at the fit (the posterior covariance of the Kalman view); the lensing style shows a sky through the recovered lens.');
      }
      ctx.setNotes(notes.join('\n'));

      function visiblePanels() {
        const sel = Math.round(options.panels);
        if (sel === 0) return panels;
        if (sel === 5) return panels.filter(p => p.name === 'observed' || p.name === 'model');
        const want = ['observed', 'initial', 'model', 'residual'][sel - 1];
        return panels.filter(p => p.name === want);
      }

      return {
        group: new THREE.Group(),
        render(r) {
          const size = r.getSize(new THREE.Vector2()), w = size.x, h = size.y;
          const pr = r.getPixelRatio();
          const rect = el.getBoundingClientRect();
          const cssScale = rect.width / w;
          const shown = visiblePanels();
          const showPlots = options.history && plots.length > 0;
          const plotH = showPlots ? Math.round(h * 0.24) : 0;
          const top = h - plotH;
          // Lay the panels out in a row, or two rows when that gives them more area.
          const n = shown.length;
          let cols = n, rows = 1;
          if (n >= 4 && (w / n) / top < 0.7) { cols = Math.ceil(n / 2); rows = 2; }
          const cellW = Math.floor(w / cols), cellH = Math.floor(top / rows);
          r.setScissorTest(true);
          r.autoClear = false;
          r.setClearColor(0x0b0d12, 1);
          r.setViewport(0, 0, w, h);
          r.setScissor(0, 0, w, h);
          r.clear();
          let k = 0;
          const ringOn = options.ring && fit && fit.fitted && fit.fitted.length >= 3;
          shown.forEach((p, idx) => {
            const cx = idx % cols, cy = Math.floor(idx / cols);
            const side = Math.min(cellW, cellH * p.aspect);
            const pw = Math.round(side - 8), ph = Math.round(side / p.aspect - 8);
            const x0 = cx * cellW + Math.round((cellW - pw) / 2), y0 = top - (cy + 1) * cellH + Math.round((cellH - ph) / 2);
            p.uniforms.uExposure.value = options.exposure;
            p.uniforms.uLog.value = options.log ? 1 : 0;
            p.uniforms.uRange.value = options.range;
            if (ringOn) p.uniforms.uRing.value.set((fit.fitted[1] - p.u0) / p.spanU, (fit.fitted[2] - p.v0) / p.spanV, fit.fitted[0]);
            else p.uniforms.uRing.value.set(0, 0, 0);
            r.setViewport(x0, y0, pw, ph);
            r.setScissor(x0, y0, pw, ph);
            r.render(p.scene, quadCamera);
            if (options.labels) {
              const d = label(k++);
              d.style.left = (rect.left + (x0 + 6) * cssScale) + 'px';
              d.style.top = (rect.top + (h - y0 - ph + 4) * cssScale) + 'px';
              d.style.display = 'block';
              d.textContent = (TITLES[p.name] || p.name) + (p.name === 'residual' ? ` (+/- ${options.range} ${p.unit || 'sigma'})` : '')
                + `\n${p.spanU.toFixed(2)}" x ${p.spanV.toFixed(2)}"`;
            }
          });
          if (showPlots) {
            r.setViewport(0, 0, w, plotH);
            r.setScissor(0, 0, w, plotH);
            r.render(plotScene, plotCamera);
            if (options.labels) plotInfo.forEach((text, i) => {
              const d = label(k++);
              d.style.left = (rect.left + (0.03 * w + 0.5 * w * i + 4) * cssScale) + 'px';
              d.style.top = (rect.top + (h - plotH + 2) * cssScale) + 'px';
              d.style.display = 'block';
              d.textContent = text;
            });
          }
          for (let i = k; i < labelEls.length; i++) labelEls[i].style.display = 'none';
          r.setScissorTest(false);
          r.setViewport(0, 0, w, h);
          r.autoClear = true;
        },
        readout() {
          if (!fit) return ['no lensFit record: images only'];
          const lines = [`Einstein radius ${fit.thetaEArcsec.toFixed(4)}"; D_l ${fit.distancesMpc[0].toFixed(0)}, D_s ${fit.distancesMpc[1].toFixed(0)}, D_ls ${fit.distancesMpc[2].toFixed(0)} Mpc`];
          const fmt = (x) => Math.abs(x) >= 1e4 || (Math.abs(x) < 1e-3 && x !== 0) ? x.toExponential(3) : x.toFixed(4);
          lines.push('parameter: first estimate -> fit +/- error' + (fit.truth && fit.truth.length ? ' (truth)' : ''));
          fit.parameterNames.forEach((name, i) => {
            const t = fit.truth && fit.truth.length > i ? ` (${fmt(fit.truth[i])})` : '';
            lines.push(`${name}: ${fmt(fit.initial[i])} -> ${fmt(fit.fitted[i])} +/- ${fmt(fit.sigma[i])}${t}`);
          });
          return lines;
        },
        onOption() { return true; },
        dispose() {
          el.removeEventListener('wheel', onWheel);
          overlay.remove();
          controls.enabled = true;
          renderer.setClearColor(previousClear, previousAlpha);
          for (const p of panels) { p.material.dispose(); p.texture.dispose(); }
          for (const l of plots) { l.geometry.dispose(); l.material.dispose(); }
        },
      };
    },
  });
})();
