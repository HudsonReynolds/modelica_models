// Lens-fit style: the results of lens inversions (GeneralRelativity.Lensing.Inversion), one case or a library of them.
// Laid out as HTML over the view, so captions, colour bars, plots and the parameter table each have their own space:
//   gallery of cases | title, description and key numbers
//                    | observed, first simulation, fitted simulation, residual (with colour bars)
//                    | cost and mass per iteration
//                    | parameters: first estimate, fit +/- error, truth
// Records: lensImage (case, name, u, v, values[i][j] at (u[i], v[j]) in arcseconds) and lensFit (case, title, description,
// parameters, errors, truth, histories), both from Visualization.writeLensInversion.
(function () {
  const ORDER = ['observed', 'initial', 'model', 'residual'];
  const TITLES = { observed: 'Observed image', initial: 'First simulation', model: 'Fitted simulation', residual: 'Residual' };
  const SUBTITLES = {
    observed: 'the input, lens light removed',
    initial: 'round lens of the ring radius',
    model: 'after the least-squares fit',
    residual: '(observed − fitted) / noise',
  };
  const CSS = `
    .lf-root { position: absolute; inset: 0; display: grid; grid-template-columns: 216px minmax(0, 1fr); background: var(--bg, #0b0e14);
      color: var(--text, #e6e9ef); font: 13px/1.4 system-ui, sans-serif; z-index: 2; outline: none; }
    .lf-root.lf-single { grid-template-columns: minmax(0, 1fr); }
    .lf-root.lf-single .lf-gallery { display: none; }
    .lf-gallery { overflow-y: auto; border-right: 1px solid var(--line, #222a38); padding: 10px 8px; display: flex; flex-direction: column; gap: 6px; min-height: 0; }
    .lf-gallery h3 { margin: 2px 4px 4px; font-size: 11px; font-weight: 600; color: var(--muted, #8a93a6); text-transform: uppercase; letter-spacing: .05em; }
    .lf-card { display: grid; grid-template-columns: 52px minmax(0, 1fr); gap: 8px; align-items: center; padding: 6px; border-radius: 6px;
      border: 1px solid transparent; background: none; color: inherit; font: inherit; text-align: left; cursor: pointer; }
    .lf-card:hover { background: #151a24; border-color: #2f3a4e; }
    .lf-card.lf-active { background: #18202e; border-color: var(--accent, #5fb3ff); }
    .lf-card canvas { width: 52px; height: 52px; border-radius: 4px; image-rendering: pixelated; background: #000; }
    .lf-card .lf-card-title { font-weight: 600; font-size: 12.5px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
    .lf-card .lf-card-sub { font-size: 11.5px; color: var(--muted, #8a93a6); }
    .lf-good { color: #7cd992; } .lf-fair { color: #f0c674; } .lf-poor { color: #ff8a7a; }
    .lf-main { overflow-y: auto; min-height: 0; padding: 14px 16px 24px; display: flex; flex-direction: column; gap: 14px; }
    .lf-head h2 { margin: 0; font-size: 17px; font-weight: 650; }
    .lf-head p { margin: 3px 0 0; color: var(--muted, #8a93a6); max-width: 80ch; }
    .lf-stats { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 10px; }
    .lf-stat { background: #151a24; border: 1px solid #222a38; border-radius: 6px; padding: 6px 10px; min-width: 0; }
    .lf-stat .lf-k { font-size: 11px; color: var(--muted, #8a93a6); }
    .lf-stat .lf-v { font-size: 14px; font-weight: 600; font-variant-numeric: tabular-nums; white-space: nowrap; }
    .lf-stat .lf-v small { font-size: 11.5px; font-weight: 500; color: var(--muted, #8a93a6); }
    .lf-panels { display: grid; gap: 12px; justify-content: start; }
    .lf-fig { margin: 0; display: flex; flex-direction: column; gap: 4px; min-width: 0; }
    .lf-fig figcaption { display: flex; align-items: baseline; gap: 8px; min-width: 0; }
    .lf-fig figcaption b { font-size: 13px; font-weight: 600; white-space: nowrap; }
    .lf-fig figcaption span { font-size: 11.5px; color: var(--muted, #8a93a6); overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
    .lf-fig .lf-img { position: relative; line-height: 0; background: #000; border-radius: 4px; overflow: hidden; }
    .lf-fig .lf-img canvas { cursor: crosshair; }
    .lf-bar { display: grid; grid-template-columns: auto 1fr auto; gap: 6px; align-items: center; font-size: 11px; color: var(--muted, #8a93a6);
      font-variant-numeric: tabular-nums; }
    .lf-bar i { display: block; height: 7px; border-radius: 3px; }
    .lf-plots { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 12px; }
    .lf-plot { background: #10141c; border: 1px solid #222a38; border-radius: 6px; padding: 8px 10px 6px; min-width: 0; }
    .lf-plot h4 { margin: 0 0 2px; font-size: 12.5px; font-weight: 600; }
    .lf-plot canvas { display: block; width: 100%; height: 170px; }
    .lf-table-wrap { overflow-x: auto; }
    .lf-table { border-collapse: collapse; font-size: 12.5px; font-variant-numeric: tabular-nums; min-width: 560px; }
    .lf-table th { text-align: right; font-weight: 600; color: var(--muted, #8a93a6); padding: 4px 10px; border-bottom: 1px solid #2f3a4e; white-space: nowrap; }
    .lf-table td { text-align: right; padding: 3px 10px; border-bottom: 1px solid #1a2030; white-space: nowrap; }
    .lf-table th:first-child, .lf-table td:first-child { text-align: left; }
    .lf-table tr.lf-group td { color: var(--muted, #8a93a6); font-size: 11px; text-transform: uppercase; letter-spacing: .05em; padding-top: 10px; border-bottom: none; }
    .lf-section-title { margin: 0; font-size: 13px; font-weight: 600; }
    .lf-tip { position: absolute; pointer-events: none; background: rgba(12, 15, 22, .94); border: 1px solid #2f3a4e; border-radius: 4px;
      padding: 3px 7px; font-size: 11.5px; font-variant-numeric: tabular-nums; white-space: nowrap; display: none; z-index: 3; }
    .lf-root.lf-narrow { grid-template-columns: minmax(0, 1fr); grid-template-rows: auto minmax(0, 1fr); }
    .lf-root.lf-narrow .lf-gallery { flex-direction: row; overflow-x: auto; overflow-y: hidden; border-right: none; border-bottom: 1px solid var(--line, #222a38); }
    .lf-root.lf-narrow .lf-gallery h3 { display: none; }
    .lf-root.lf-narrow .lf-card { grid-template-columns: 40px 120px; flex: none; }
    .lf-root.lf-narrow .lf-card canvas { width: 40px; height: 40px; }
    /* The viewer stretches every canvas in #view over the whole view (for the WebGL canvas); undo that inside this layout.
       These selectors include the id so they outrank that rule. */
    #view .lf-root canvas { position: static; inset: auto; width: auto; height: auto; display: block; }
    #view .lf-root .lf-card canvas { width: 52px; height: 52px; }
    #view .lf-root.lf-narrow .lf-card canvas { width: 40px; height: 40px; }
    #view .lf-root .lf-plot canvas { width: 100%; height: 170px; }
    #view .lf-root .lf-img canvas { position: absolute; inset: 0; width: 100%; height: 100%; }
  `;

  // Colour maps: brightness (black - purple - orange - cream) and residual (blue - grey - red).
  const HEAT = [[5, 3, 20], [115, 13, 115], [250, 115, 26], [255, 255, 217]];
  function heat(t) {
    t = Math.min(Math.max(t, 0), 1) * 3;
    const k = Math.min(Math.floor(t), 2), f = t - k, a = HEAT[k], b = HEAT[k + 1];
    return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f];
  }
  function diverging(t) {
    t = Math.min(Math.max(t, -1), 1);
    const mid = [34, 36, 44], end = t < 0 ? [64, 140, 255] : [255, 90, 51], f = Math.abs(t);
    return [mid[0] + (end[0] - mid[0]) * f, mid[1] + (end[1] - mid[1]) * f, mid[2] + (end[2] - mid[2]) * f];
  }
  const rgb = (c) => `rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0})`;
  const gradient = (fn, a, b) => 'linear-gradient(90deg,' + Array.from({ length: 9 }, (_, k) => rgb(fn(a + (b - a) * k / 8))).join(',') + ')';
  const SUP = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
  const sup = (n) => String(n).split('').map(c => SUP[c] || c).join('');
  // A value with its error in one power of ten, with as many digits as the error warrants: 2.872 +/- 0.002 x 10^11.
  function withError(value, sigma) {
    const e = Math.floor(Math.log10(Math.abs(value) || 1)), scale = Math.pow(10, e);
    const d = sigma > 0 ? Math.min(Math.max(-Math.floor(Math.log10(sigma / scale)), 1), 4) : 3;
    return { text: (value / scale).toFixed(d) + (sigma > 0 ? ' ± ' + (sigma / scale).toFixed(d) : ''), unit: '× 10' + sup(e) };
  }
  const sci = (x, d = 3) => {
    if (!isFinite(x)) return '—';
    if (x === 0) return '0';
    const e = Math.floor(Math.log10(Math.abs(x)));
    return Math.abs(e) < 4 ? x.toFixed(Math.max(0, d - 1 - e)) : (x / Math.pow(10, e)).toFixed(d - 1) + ' × 10' + sup(e);
  };
  const num = (x, d = 4) => (!isFinite(x) ? '—' : Math.abs(x) >= 1e4 || (Math.abs(x) < 1e-3 && x !== 0) ? sci(x, 3) : x.toFixed(d));
  function h(tag, attrs, ...children) {
    const el = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs || {})) {
      if (k === 'class') el.className = v; else if (k === 'style') el.style.cssText = v; else if (k.startsWith('on')) el.addEventListener(k.slice(2), v); else el.setAttribute(k, v);
    }
    for (const c of children.flat()) if (c != null) el.append(c instanceof Node ? c : document.createTextNode(String(c)));
    return el;
  }

  // Cases: every lensFit record, then any images whose case has no fit record.
  function casesOf(data) {
    const recs = data.records || {}, fits = recs.lensFit || [], images = recs.lensImage || [];
    const cases = [], byKey = new Map();
    const get = (key) => {
      if (!byKey.has(key)) { const c = { key, fit: null, images: {} }; byKey.set(key, c); cases.push(c); }
      return byKey.get(key);
    };
    for (const f of fits) get(f.case || '').fit = f;
    for (const im of images) get(im.case || '').images[im.name] = im;
    for (const c of cases) c.title = (c.fit && c.fit.title) || c.key || (data.header && data.header.name) || 'Lens inversion';
    return cases;
  }
  // How well the fit recovered the truth, for the gallery and the header.
  function grade(fit) {
    if (!fit) return { text: 'images only', cls: '' };
    const last = fit.chi2History[Math.min(fit.iterations, fit.chi2History.length - 1)];
    const chi = fit.pixels > 0 ? last / fit.pixels : NaN;
    if (fit.trueMassSolar > 0 && fit.massSigmaSolar > 0) {
      const z = (fit.massSolar - fit.trueMassSolar) / fit.massSigmaSolar, a = Math.abs(z);
      return { text: `${z >= 0 ? '+' : '−'}${a.toFixed(1)}σ from truth`, cls: a < 2 ? 'lf-good' : a < 4 ? 'lf-fair' : 'lf-poor', z, chi };
    }
    return { text: `χ²/N ${chi.toFixed(2)}`, cls: chi < 1.3 ? 'lf-good' : chi < 2 ? 'lf-fair' : 'lf-poor', chi };
  }

  // A small line plot on a 2D canvas: axes with nice ticks, optional log y, dashed reference lines, a shaded band and a legend.
  function niceTicks(lo, hi, n) {
    const span = hi - lo, raw = span / Math.max(n, 1), mag = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / mag;
    const step = (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * mag, out = [];
    for (let t = Math.ceil(lo / step) * step; t <= hi + 1e-9 * step; t += step) out.push(t);
    return out;
  }
  function drawPlot(canvas, cfg, colors) {
    const dpr = window.devicePixelRatio || 1, W = canvas.clientWidth, H = canvas.clientHeight;
    if (!W || !H) return;
    canvas.width = Math.round(W * dpr); canvas.height = Math.round(H * dpr);
    const g = canvas.getContext('2d');
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.clearRect(0, 0, W, H);
    const m = { l: 54, r: 10, t: 22, b: 30 }, pw = W - m.l - m.r, ph = H - m.t - m.b;
    const ys = cfg.series.flatMap(s => s.points.map(p => p[1])).concat(cfg.lines.map(l => l.y)).concat(cfg.band ? cfg.band : []).filter(v => isFinite(v) && (!cfg.log || v > 0));
    let y0 = Math.min(...ys), y1 = Math.max(...ys);
    if (cfg.log) { y0 = Math.log10(y0); y1 = Math.log10(y1); }
    const pad = Math.max((y1 - y0) * 0.08, cfg.log ? 0.05 : Math.abs(y1) * 1e-3 || 1e-3);
    y0 -= pad; y1 += pad;
    const xmax = Math.max(1, cfg.xmax);
    const X = (x) => m.l + pw * x / xmax, Y = (y) => m.t + ph * (1 - ((cfg.log ? Math.log10(y) : y) - y0) / (y1 - y0));
    g.font = '11px system-ui, sans-serif';
    g.strokeStyle = colors.grid; g.fillStyle = colors.muted; g.lineWidth = 1;
    // y ticks
    let yt;
    if (cfg.log) {
      yt = [];
      for (let e = Math.floor(y0); e <= Math.ceil(y1); e++) for (const k of (y1 - y0 < 1.5 ? [1, 2, 5] : [1])) { const v = k * Math.pow(10, e); const lv = Math.log10(v); if (lv >= y0 && lv <= y1) yt.push(v); }
    } else yt = niceTicks(y0, y1, 4);
    g.textAlign = 'right'; g.textBaseline = 'middle';
    for (const v of yt) {
      const y = Y(v);
      g.beginPath(); g.moveTo(m.l, y + 0.5); g.lineTo(m.l + pw, y + 0.5); g.stroke();
      g.fillText(cfg.yfmt(v), m.l - 6, y);
    }
    g.textAlign = 'center'; g.textBaseline = 'top';
    for (const v of niceTicks(0, xmax, 6)) {
      const x = X(v);
      g.beginPath(); g.moveTo(x + 0.5, m.t + ph); g.lineTo(x + 0.5, m.t + ph + 4); g.stroke();
      g.fillText(String(Math.round(v)), x, m.t + ph + 6);
    }
    g.strokeStyle = colors.axis;
    g.beginPath(); g.moveTo(m.l + 0.5, m.t); g.lineTo(m.l + 0.5, m.t + ph + 0.5); g.lineTo(m.l + pw, m.t + ph + 0.5); g.stroke();
    g.fillStyle = colors.muted; g.textAlign = 'right'; g.textBaseline = 'bottom';
    g.fillText(cfg.xlabel, m.l + pw, H - 1);
    g.textAlign = 'left'; g.textBaseline = 'top';
    g.fillText(cfg.ylabel, 4, 3);
    if (cfg.band) {
      g.fillStyle = cfg.bandColor;
      const ya = Y(cfg.band[0]), yb = Y(cfg.band[1]);
      g.fillRect(m.l, Math.min(ya, yb), pw, Math.max(Math.abs(yb - ya), 1));
    }
    g.save();
    g.beginPath(); g.rect(m.l, m.t - 2, pw, ph + 4); g.clip();
    for (const l of cfg.lines) {
      g.strokeStyle = l.color; g.setLineDash(l.dash || [5, 4]); g.lineWidth = 1.5;
      g.beginPath(); g.moveTo(m.l, Y(l.y)); g.lineTo(m.l + pw, Y(l.y)); g.stroke();
    }
    g.setLineDash([]);
    for (const s of cfg.series) {
      g.strokeStyle = s.color; g.fillStyle = s.color; g.lineWidth = 2;
      if (s.points.length > 1) {
        g.beginPath();
        s.points.forEach((p, i) => (i ? g.lineTo(X(p[0]), Y(p[1])) : g.moveTo(X(p[0]), Y(p[1]))));
        g.stroke();
      }
      for (const p of s.markers ? s.points : []) { g.beginPath(); g.arc(X(p[0]), Y(p[1]), s.radius || 2.2, 0, 2 * Math.PI); g.fill(); }
    }
    g.restore();
    // Legend, top right, inside the margin row.
    const items = cfg.series.filter(s => s.label).map(s => ({ label: s.label, color: s.color, dash: null }))
      .concat(cfg.lines.filter(l => l.label).map(l => ({ label: l.label, color: l.color, dash: l.dash || [5, 4] })))
      .concat(cfg.band && cfg.bandLabel ? [{ label: cfg.bandLabel, color: cfg.bandColor, box: true }] : []);
    let x = W - m.r;
    g.textAlign = 'right'; g.textBaseline = 'middle';
    for (const it of items.reverse()) {
      g.fillStyle = colors.text;
      g.fillText(it.label, x, 9);
      x -= g.measureText(it.label).width + 6;
      if (it.box) { g.fillStyle = it.color; g.fillRect(x - 14, 5, 14, 8); }
      else { g.strokeStyle = it.color; g.lineWidth = 2; g.setLineDash(it.dash || []); g.beginPath(); g.moveTo(x - 16, 9); g.lineTo(x, 9); g.stroke(); g.setLineDash([]); }
      x -= 28;
    }
  }

  window.GRStyles = window.GRStyles || [];
  window.GRStyles.push({
    id: 'lensfit',
    label: 'Lens inversion',
    order: 8,
    description: 'The mass behind lensed images: pick a case in the gallery to see the observed image, the simulation of the first estimate from the ring radius, the fitted simulation and the residual, with the cost and mass per iteration and every fitted parameter. Hover an image for pixel values; arrow keys step through the cases.',
    options: [
      { id: 'case', label: 'case', type: 'select', value: 0, choices: (data) => casesOf(data).map((c, i) => ({ value: i, label: c.title })) },
      { id: 'panels', label: 'panels', type: 'select', choices: [{ value: 0, label: 'all four' }, { value: 1, label: 'observed and fitted' }, { value: 2, label: 'observed, fitted, residual' }], value: 0 },
      { id: 'exposure', label: 'exposure', type: 'range', min: 0.2, max: 8, step: 0.1, value: 1.2 },
      { id: 'log', label: 'log stretch', type: 'checkbox', value: 1 },
      { id: 'range', label: 'residual range (σ)', type: 'range', min: 1, max: 20, step: 0.5, value: 4 },
      { id: 'ring', label: 'critical curve', type: 'checkbox', value: 1, help: 'The fitted lens\'s critical curve (dashed) and centre (+); the true centre (×) when known' },
      { id: 'plots', label: 'iteration plots', type: 'checkbox', value: 1 },
      { id: 'table', label: 'parameter table', type: 'checkbox', value: 1 },
    ],
    available(data) {
      return data.records && data.records.lensImage ? null : 'needs lens-inversion images (Examples.LensInversion, LensInversionLibrary or Visualization.writeLensInversion)';
    },
    create(ctx) {
      const { THREE, data, options, renderer, controls } = ctx;
      controls.enabled = false;
      const previousClear = renderer.getClearColor(new THREE.Color()), previousAlpha = renderer.getClearAlpha();
      const cases = casesOf(data);
      options.case = Math.min(Math.max(Math.round(options.case || 0), 0), cases.length - 1);

      const host = renderer.domElement.parentElement;
      let styleTag = document.getElementById('lensfit-css');
      if (!styleTag) { styleTag = h('style', { id: 'lensfit-css' }, CSS); document.head.appendChild(styleTag); }
      const css = getComputedStyle(document.documentElement);
      const colors = {
        text: css.getPropertyValue('--text').trim() || '#e6e9ef', muted: css.getPropertyValue('--muted').trim() || '#8a93a6',
        grid: '#1f2633', axis: '#3a4458', accent: css.getPropertyValue('--accent').trim() || '#5fb3ff',
      };

      const gallery = h('nav', { class: 'lf-gallery', 'aria-label': 'Cases' }, h('h3', null, `${cases.length} cases`));
      const main = h('div', { class: 'lf-main' });
      const tip = h('div', { class: 'lf-tip' });
      const root = h('div', { class: 'lf-root' + (cases.length > 1 ? '' : ' lf-single'), tabindex: '0' }, gallery, main, tip);
      host.appendChild(root);

      // Gallery cards with thumbnails of the observed images.
      const cards = cases.map((c, i) => {
        const thumb = h('canvas', { width: '64', height: '64' });
        const obs = c.images.observed || Object.values(c.images)[0];
        if (obs) paint(thumb, obs, 'brightness', imagePeak(obs), 1.4, true, 4);
        const gr = grade(c.fit);
        const card = h('button', { class: 'lf-card', type: 'button', title: (c.fit && c.fit.description) || c.title },
          thumb, h('div', null, h('div', { class: 'lf-card-title' }, c.title), h('div', { class: 'lf-card-sub ' + gr.cls }, gr.text)));
        card.addEventListener('click', () => selectCase(i));
        gallery.appendChild(card);
        return card;
      });
      root.addEventListener('keydown', (e) => {
        if (e.key === 'ArrowDown' || e.key === 'ArrowRight') { selectCase(Math.min(options.case + 1, cases.length - 1)); e.preventDefault(); }
        if (e.key === 'ArrowUp' || e.key === 'ArrowLeft') { selectCase(Math.max(options.case - 1, 0)); e.preventDefault(); }
      });
      function selectCase(i) {
        if (i === options.case) return;
        options.case = i;
        build();
        ctx.refreshOptionsUI();
        ctx.refreshReadout();
        cards[i].scrollIntoView({ block: 'nearest', inline: 'nearest' });
      }

      function imagePeak(rec) {
        let p = 0;
        for (const row of rec.values) for (const v of row) if (v > p) p = v;
        return p || 1;
      }
      // Draw an image record into a canvas at its native resolution (CSS scales it up with pixelated rendering).
      function paint(canvas, rec, mode, peak, exposure, log, range) {
        const nu = rec.u.length, nv = rec.v.length;
        canvas.width = nu; canvas.height = nv;
        const g = canvas.getContext('2d'), img = g.createImageData(nu, nv);
        for (let i = 0; i < nu; i++) for (let j = 0; j < nv; j++) {
          const val = rec.values[i][j];
          let c;
          if (mode === 'residual') c = diverging(val / range);
          else { const s = Math.max(val, 0) / peak * exposure; c = heat(log ? Math.log1p(30 * s) / Math.log(31) : s); }
          const k = ((nv - 1 - j) * nu + i) * 4;
          img.data[k] = c[0]; img.data[k + 1] = c[1]; img.data[k + 2] = c[2]; img.data[k + 3] = 255;
        }
        g.putImageData(img, 0, 0);
      }

      let panels = [], plotCanvases = [], tile = 200, observer = null;
      function build() {
        cards.forEach((c, i) => c.classList.toggle('lf-active', i === options.case));
        main.textContent = '';
        const c = cases[options.case], fit = c.fit;
        // Header: title, description, key numbers.
        const head = h('div', { class: 'lf-head' }, h('h2', null, c.title), fit && fit.description ? h('p', null, fit.description) : null);
        if (fit) {
          const gr = grade(fit);
          const m = withError(fit.massSolar, fit.massSigmaSolar);
          const stat = (k, v, small, cls) => h('div', { class: 'lf-stat' }, h('div', { class: 'lf-k' }, k), h('div', { class: 'lf-v ' + (cls || '') }, v, small ? h('small', null, ' ' + small) : null));
          const stats = [
            stat('Mass inside the Einstein radius', m.text, m.unit + ' M☉'),
            stat('Estimate from the ring radius', sci(fit.massInitialSolar), 'M☉'),
          ];
          if (fit.trueMassSolar > 0) {
            const pct = 100 * (fit.massSolar / fit.trueMassSolar - 1);
            stats.push(stat('True mass', sci(fit.trueMassSolar), `M☉ · fit ${pct >= 0 ? '+' : '−'}${Math.abs(pct).toFixed(2)}%`));
            stats.push(stat('Accuracy', gr.text, '', gr.cls));
          }
          stats.push(stat('Einstein radius', fit.thetaEArcsec.toFixed(3) + '″', `lens at ${fit.distancesMpc[0].toFixed(0)} Mpc`));
          stats.push(stat('Fit quality', `χ²/N ${gr.chi.toFixed(3)}`, `${fit.iterations} iterations${fit.converged ? '' : ', limit reached'}`, gr.chi < 1.3 ? 'lf-good' : gr.chi < 2 ? 'lf-fair' : 'lf-poor'));
          stats.push(stat('Image explained', (100 * fit.explained).toFixed(1) + '%', `correlation ${fit.correlation.toFixed(3)}`));
          head.appendChild(h('div', { class: 'lf-stats' }, stats));
        }
        main.appendChild(head);

        // Image panels.
        const sel = Math.round(options.panels);
        const wanted = sel === 1 ? ['observed', 'model'] : sel === 2 ? ['observed', 'model', 'residual'] : ORDER;
        const names = wanted.filter(n => c.images[n]).concat(sel === 0 ? Object.keys(c.images).filter(n => !ORDER.includes(n)) : []);
        const peak = c.images.observed ? imagePeak(c.images.observed) : Math.max(...names.filter(n => n !== 'residual').map(n => imagePeak(c.images[n])), 1e-30);
        const grid = h('section', { class: 'lf-panels' });
        panels = names.map(name => {
          const rec = c.images[name], residual = name === 'residual' || rec.unit === 'sigma';
          const nu = rec.u.length, nv = rec.v.length, du = nu > 1 ? rec.u[1] - rec.u[0] : 1, dv = nv > 1 ? rec.v[1] - rec.v[0] : 1;
          const geom = { u0: rec.u[0] - du / 2, v0: rec.v[0] - dv / 2, spanU: nu * du, spanV: nv * dv };
          const base = h('canvas', { style: 'image-rendering: pixelated; position: absolute; inset: 0; width: 100%; height: 100%;' });
          const over = h('canvas', { style: 'position: absolute; inset: 0; width: 100%; height: 100%;' });
          const box = h('div', { class: 'lf-img' }, base, over);
          paint(base, rec, residual ? 'residual' : 'brightness', peak, options.exposure, options.log, options.range);
          const bar = residual
            ? h('div', { class: 'lf-bar' }, `−${options.range}σ`, h('i', { style: 'background:' + gradient(diverging, -1, 1) }), `+${options.range}σ`)
            : h('div', { class: 'lf-bar' }, '0', h('i', { style: 'background:' + gradient(t => heat(t), 0, 1) }), `${sci(peak / options.exposure, 2)}${options.log ? ' (log)' : ''}`);
          const fig = h('figure', { class: 'lf-fig' },
            h('figcaption', null, h('b', null, TITLES[name] || name), h('span', null, SUBTITLES[name] || '')), box, bar);
          grid.appendChild(fig);
          const p = { name, rec, residual, geom, base, over, box, aspect: geom.spanU / geom.spanV };
          over.addEventListener('mousemove', (e) => hover(p, e));
          over.addEventListener('mouseleave', () => { tip.style.display = 'none'; });
          return p;
        });
        main.appendChild(grid);

        // Plots.
        plotCanvases = [];
        if (fit && options.plots && fit.chi2History && fit.chi2History.length > 1) {
          const n = Math.min(fit.iterations, fit.chi2History.length - 1);
          const N = fit.pixels > 0 ? fit.pixels : 1;
          const chiCanvas = h('canvas'), massCanvas = h('canvas');
          main.appendChild(h('section', { class: 'lf-plots' },
            h('div', { class: 'lf-plot' }, h('h4', null, 'Fit cost per iteration'), chiCanvas),
            h('div', { class: 'lf-plot' }, h('h4', null, 'Mass estimate per iteration'), massCanvas)));
          const chi = fit.chi2History.slice(0, n + 1).map((v, i) => [i, v / N]);
          plotCanvases.push([chiCanvas, {
            xmax: n, log: true, xlabel: 'iteration', ylabel: 'χ² / pixels', yfmt: (v) => (v >= 100 ? v.toExponential(0).replace('e+', 'e') : String(+v.toPrecision(2))),
            series: [{ points: chi, color: '#ffb454', markers: n <= 60, label: 'simulation vs image' }],
            lines: [{ y: 1, color: '#8a93a6', label: 'noise level' }],
          }]);
          const mh = (fit.massHistorySolar || []).slice(0, n + 1);
          if (mh.length) {
            const all = mh.concat(fit.trueMassSolar > 0 ? [fit.trueMassSolar] : []).concat([fit.massInitialSolar]);
            const e = Math.floor(Math.log10(Math.max(...all))), s = Math.pow(10, e);
            plotCanvases.push([massCanvas, {
              xmax: n, log: false, xlabel: 'iteration', ylabel: `mass (10${sup(e)} M☉)`, yfmt: (v) => String(+v.toPrecision(3)),
              series: [
                { points: mh.map((m, i) => [i, m / s]), color: '#6ad0ff', markers: n <= 60, label: 'fit' },
                { points: [[0, fit.massInitialSolar / s]], color: '#ffb454', markers: true, radius: 4, label: 'ring radius' },
              ],
              lines: fit.trueMassSolar > 0 ? [{ y: fit.trueMassSolar / s, color: '#7cd992', label: 'truth' }] : [],
              band: fit.massSigmaSolar > 0 ? [(fit.massSolar - fit.massSigmaSolar) / s, (fit.massSolar + fit.massSigmaSolar) / s] : null,
              bandColor: 'rgba(106, 208, 255, 0.22)', bandLabel: '±1σ',
            }]);
          }
        }

        // Parameter table.
        if (fit && options.table && fit.parameterNames) {
          const truth = fit.truth && fit.truth.length === fit.parameterNames.length ? fit.truth : null;
          const rows = [];
          fit.parameterNames.forEach((name, i) => {
            if (i === 0) rows.push(h('tr', { class: 'lf-group' }, h('td', { colspan: truth ? '6' : '4' }, 'Lens')));
            if (i === 7) rows.push(h('tr', { class: 'lf-group' }, h('td', { colspan: truth ? '6' : '4' }, 'Source')));
            const sig = fit.sigma[i], fixed = !(sig > 0);
            let off = null;
            if (truth) {
              const z = (fit.fitted[i] - truth[i]) / sig;
              off = fixed ? h('td', { class: 'lf-k' }, 'fixed') : h('td', { class: Math.abs(z) < 2 ? 'lf-good' : Math.abs(z) < 4 ? 'lf-fair' : 'lf-poor' }, `${z >= 0 ? '+' : '−'}${Math.abs(z).toFixed(1)}σ`);
            }
            rows.push(h('tr', null, h('td', null, name), h('td', null, num(fit.initial[i])), h('td', null, num(fit.fitted[i])),
              h('td', null, fixed ? '—' : '± ' + num(sig)), truth ? h('td', null, num(truth[i])) : null, off));
          });
          const header = h('tr', null, h('th', null, 'Parameter'), h('th', null, 'First estimate'), h('th', null, 'Fit'), h('th', null, 'Error'),
            truth ? h('th', null, 'Truth') : null, truth ? h('th', null, 'Fit − truth') : null);
          main.appendChild(h('section', null, h('h3', { class: 'lf-section-title' }, 'Parameters'),
            h('div', { class: 'lf-table-wrap' }, h('table', { class: 'lf-table' }, h('thead', null, header), h('tbody', null, rows)))));
        }

        const notes = [];
        if (fit) {
          notes.push('Mass is the projected mass inside the critical curve, πq b² Σcr Dₗ² (Lensing.Inversion.einsteinMass); errors come from the curvature of χ² at the fit, the posterior covariance of the Kalman view.');
          if (!fit.converged) notes.push('This fit stopped at the iteration limit; raise maxIterations for a fully converged answer.');
        }
        ctx.setNotes(notes.join('\n'));
        layout();
      }

      // Size the image tiles so that the header, the images and the plots fit the view height where possible.
      function layout() {
        const W = main.clientWidth - 32;
        if (!W || !panels.length) return;
        root.classList.toggle('lf-narrow', root.clientWidth < 720 && cases.length > 1);
        const head = main.querySelector('.lf-head'), plots = main.querySelector('.lf-plots');
        const reserved = (head ? head.offsetHeight : 0) + (plots ? plots.offsetHeight + 14 : 0) + 14 + 38 + 40;
        const Hh = Math.max(main.clientHeight - reserved, 160);
        const n = panels.length, gap = 12, aspect = panels[0].aspect;
        // The largest tile that fits both the width and the height. When that is too small to read, give up on fitting the height:
        // use the fewest columns whose tiles are at most 260 px wide, and let the page scroll.
        const byWidth = (cols) => (W - gap * (cols - 1)) / cols;
        let best = 0, bestCols = n;
        for (let cols = 1; cols <= n; cols++) {
          const rows = Math.ceil(n / cols);
          const s = Math.min(byWidth(cols), ((Hh - gap * (rows - 1)) / rows - 44) * aspect);
          if (s > best + 0.5) { best = s; bestCols = cols; }
        }
        if (best < 160) {
          bestCols = n;
          for (let cols = 1; cols <= n; cols++) if (byWidth(cols) <= 260) { bestCols = cols; break; }
          best = Math.min(byWidth(bestCols), 260);
        }
        tile = Math.max(Math.floor(Math.min(best, 520)), 60);
        const grid = main.querySelector('.lf-panels');
        grid.style.gridTemplateColumns = `repeat(${bestCols}, ${tile}px)`;
        for (const p of panels) {
          p.box.style.width = tile + 'px';
          p.box.style.height = Math.round(tile / p.aspect) + 'px';
          drawOverlay(p);
        }
        for (const [canvas, cfg] of plotCanvases) drawPlot(canvas, cfg, colors);
      }

      // Critical curve and centre of the fitted lens (and the true centre) over each image.
      function drawOverlay(p) {
        const c = cases[options.case], fit = c.fit, cv = p.over, dpr = window.devicePixelRatio || 1;
        const W = p.box.clientWidth, H = p.box.clientHeight;
        cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr);
        const g = cv.getContext('2d');
        g.setTransform(dpr, 0, 0, dpr, 0, 0);
        g.clearRect(0, 0, W, H);
        if (!fit || !options.ring || !fit.fitted) return;
        const X = (x) => (x - p.geom.u0) / p.geom.spanU * W, Y = (y) => (1 - (y - p.geom.v0) / p.geom.spanV) * H;
        const [b, lx, ly, q, paDeg] = fit.fitted, point = fit.model === 'point', pa = paDeg * Math.PI / 180;
        g.strokeStyle = p.residual ? 'rgba(255,255,255,0.75)' : 'rgba(110, 230, 255, 0.85)';
        g.lineWidth = 1.25; g.setLineDash([5, 4]);
        g.beginPath();
        // Critical curve of the isothermal ellipse: semi-axes b (major axis at the position angle) and q b; a circle for a point mass.
        const a1 = b, a2 = point ? b : q * b;
        for (let k = 0; k <= 96; k++) {
          const t = 2 * Math.PI * k / 96, ex = a1 * Math.cos(t), ey = a2 * Math.sin(t);
          const x = lx + Math.cos(pa) * ex - Math.sin(pa) * ey, y = ly + Math.sin(pa) * ex + Math.cos(pa) * ey;
          if (k) g.lineTo(X(x), Y(y)); else g.moveTo(X(x), Y(y));
        }
        g.stroke(); g.setLineDash([]);
        const cross = (x, y, r, diag, color) => {
          g.strokeStyle = color; g.lineWidth = 1.5; g.beginPath();
          if (diag) { g.moveTo(x - r, y - r); g.lineTo(x + r, y + r); g.moveTo(x - r, y + r); g.lineTo(x + r, y - r); }
          else { g.moveTo(x - r, y); g.lineTo(x + r, y); g.moveTo(x, y - r); g.lineTo(x, y + r); }
          g.stroke();
        };
        cross(X(lx), Y(ly), 6, false, 'rgba(110, 230, 255, 0.95)');
        if (fit.truth && fit.truth.length >= 3) cross(X(fit.truth[1]), Y(fit.truth[2]), 4.5, true, 'rgba(124, 217, 146, 0.95)');
      }

      function hover(p, e) {
        const r = p.over.getBoundingClientRect(), rr = root.getBoundingClientRect();
        const fx = (e.clientX - r.left) / r.width, fy = (e.clientY - r.top) / r.height;
        const nu = p.rec.u.length, nv = p.rec.v.length;
        const i = Math.min(Math.max(Math.floor(fx * nu), 0), nu - 1), j = Math.min(Math.max(Math.floor((1 - fy) * nv), 0), nv - 1);
        const val = p.rec.values[i][j];
        tip.textContent = `x ${p.rec.u[i].toFixed(2)}″  y ${p.rec.v[j].toFixed(2)}″   ${p.residual ? (val >= 0 ? '+' : '−') + Math.abs(val).toFixed(2) + 'σ' : sci(val, 3)}`;
        tip.style.display = 'block';
        const tx = e.clientX - rr.left + 14, ty = e.clientY - rr.top + 14;
        tip.style.left = Math.min(tx, rr.width - tip.offsetWidth - 6) + 'px';
        tip.style.top = Math.min(ty, rr.height - tip.offsetHeight - 6) + 'px';
      }

      build();
      observer = new ResizeObserver(() => layout());
      observer.observe(root);
      setTimeout(() => root.focus({ preventScroll: true }), 0);

      return {
        group: new THREE.Group(),
        render(r) {
          // The page is HTML over the view; keep the WebGL canvas plain.
          r.setClearColor(0x0b0e14, 1);
          r.clear();
        },
        readout() {
          const c = cases[options.case], fit = c && c.fit;
          if (!fit) return ['images only (no lensFit record)'];
          const m = withError(fit.massSolar, fit.massSigmaSolar), gr = grade(fit);
          const lines = [`${c.title}`, `mass ${m.text} ${m.unit} M_sun`, `ring estimate ${sci(fit.massInitialSolar)} M_sun`];
          if (fit.trueMassSolar > 0) lines.push(`truth ${sci(fit.trueMassSolar)} M_sun (${gr.text})`);
          lines.push(`theta_E ${fit.thetaEArcsec.toFixed(4)}", chi2/N ${gr.chi.toFixed(3)}`);
          return lines;
        },
        onOption(id) {
          if (id === 'case') { options.case = Math.min(Math.max(Math.round(options.case), 0), cases.length - 1); build(); return true; }
          build();
          return true;
        },
        dispose() {
          if (observer) observer.disconnect();
          root.remove();
          controls.enabled = true;
          renderer.setClearColor(previousClear, previousAlpha);
        },
      };
    },
  });
})();
