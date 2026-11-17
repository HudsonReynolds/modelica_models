// Shared helper for styles that draw 2D plots (the cosmology diagrams and the CMB views). The plot is drawn
// with the Canvas 2D API on a canvas the size of the renderer's drawing buffer, then shown through a
// full-screen textured quad, so it is part of the WebGL image (and of offline renders).
window.GRPlot = {
  // A drawing layer: begin(renderer) sizes the canvas and returns {w, h, u} (u: a length unit for fonts
  // and lines, about 1 per 600 device pixels of the smaller side); render(renderer, over) shows it, over the
  // current frame if `over`.
  layer(ctx) {
    const { THREE } = ctx;
    const canvas = document.createElement('canvas');
    canvas.width = canvas.height = 2;
    const texture = new THREE.CanvasTexture(canvas);
    texture.colorSpace = THREE.SRGBColorSpace;
    texture.generateMipmaps = false;
    texture.minFilter = texture.magFilter = THREE.LinearFilter;
    const material = new THREE.MeshBasicMaterial({ map: texture, transparent: true, depthTest: false, depthWrite: false, toneMapped: false });
    const scene = new THREE.Scene();
    scene.add(new THREE.Mesh(new THREE.PlaneGeometry(2, 2), material));
    const camera = new THREE.OrthographicCamera(-1, 1, 1, -1, 0, 1);
    const size = new THREE.Vector2();
    return {
      canvas,
      g: canvas.getContext('2d'),
      begin(r) {
        r.getDrawingBufferSize(size);
        if (canvas.width !== size.x || canvas.height !== size.y) { canvas.width = size.x; canvas.height = size.y; }
        return { w: canvas.width, h: canvas.height, u: Math.max(0.6, Math.min(canvas.width, canvas.height) / 600) };
      },
      render(r, over = false) {
        texture.needsUpdate = true;
        const autoClear = r.autoClear;
        r.autoClear = !over;
        r.setRenderTarget(null);
        r.render(scene, camera);
        r.autoClear = autoClear;
      },
      dispose() { texture.dispose(); material.dispose(); },
    };
  },

  // About n "nice" tick values in [lo, hi].
  ticks(lo, hi, n = 6) {
    const span = hi - lo;
    if (!(span > 0)) return [lo];
    const raw = span / n, mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const step = [1, 2, 2.5, 5, 10].map((s) => s * mag).find((s) => s >= raw) || 10 * mag;
    const out = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9 * span; v += step) out.push(Math.abs(v) < 1e-12 * span ? 0 : v);
    return out;
  },
  logTicks(lo, hi) {
    const out = [];
    for (let e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(hi)); e++) { const v = Math.pow(10, e); if (v >= lo * 0.999 && v <= hi * 1.001) out.push(v); }
    return out;
  },
  label(v) {
    const a = Math.abs(v);
    if (a === 0) return '0';
    if (a >= 1e5 || a < 1e-3) { const e = Math.floor(Math.log10(a)); const m = v / Math.pow(10, e); return (Math.abs(m - 1) < 1e-9 ? '' : +m.toPrecision(2) + '×') + '10' + String(e).replace(/-/g, '⁻').replace(/\d/g, (d) => '⁰¹²³⁴⁵⁶⁷⁸⁹'[d]); }
    return String(+v.toPrecision(4));
  },

  // Axes in the box {x, y, w, h} (device pixels). ax = {lo, hi, log?, map?, ticks?, label}; map(v) -> 0..1
  // overrides the scale. Returns {X(v), Y(v)} mapping data to pixels. Grid lines are faint.
  axes(g, box, ax, ay, u, colours = {}) {
    const fx = ax.map || (ax.log ? (v) => Math.log(v / ax.lo) / Math.log(ax.hi / ax.lo) : (v) => (v - ax.lo) / (ax.hi - ax.lo));
    const fy = ay.map || (ay.log ? (v) => Math.log(v / ay.lo) / Math.log(ay.hi / ay.lo) : (v) => (v - ay.lo) / (ay.hi - ay.lo));
    const X = (v) => box.x + box.w * fx(v), Y = (v) => box.y + box.h * (1 - fy(v));
    const text = colours.text || '#c9d1e0', grid = colours.grid || 'rgba(255,255,255,0.08)';
    g.save();
    g.font = `${Math.round(11 * u)}px system-ui, sans-serif`;
    g.lineWidth = Math.max(1, u);
    g.strokeStyle = grid;
    g.fillStyle = text;
    const xt = ax.ticks || (ax.log ? this.logTicks(ax.lo, ax.hi) : this.ticks(ax.lo, ax.hi, Math.max(3, Math.round(box.w / (90 * u)))));
    const yt = ay.ticks || (ay.log ? this.logTicks(ay.lo, ay.hi) : this.ticks(ay.lo, ay.hi, Math.max(3, Math.round(box.h / (55 * u)))));
    g.textAlign = 'center'; g.textBaseline = 'top';
    for (const v of xt) { const x = X(v); g.beginPath(); g.moveTo(x, box.y); g.lineTo(x, box.y + box.h); g.stroke(); if (!ax.noLabels) g.fillText((ax.format || this.label)(v), x, box.y + box.h + 4 * u); }
    g.textAlign = 'right'; g.textBaseline = 'middle';
    for (const v of yt) { const y = Y(v); g.beginPath(); g.moveTo(box.x, y); g.lineTo(box.x + box.w, y); g.stroke(); g.fillText((ay.format || this.label)(v), box.x - 5 * u, y); }
    g.strokeStyle = 'rgba(255,255,255,0.35)';
    g.strokeRect(box.x, box.y, box.w, box.h);
    g.font = `${Math.round(12 * u)}px system-ui, sans-serif`;
    if (ax.label) { g.textAlign = 'center'; g.textBaseline = 'top'; g.fillText(ax.label, box.x + box.w / 2, box.y + box.h + 20 * u); }
    if (ay.label) { g.save(); g.translate(box.x - 44 * u, box.y + box.h / 2); g.rotate(-Math.PI / 2); g.textAlign = 'center'; g.textBaseline = 'bottom'; g.fillText(ay.label, 0, 0); g.restore(); }
    g.restore();
    return { X, Y };
  },

  // Polyline through (xs[i], ys[i]) with the mapping {X, Y}, clipped to the box.
  line(g, map, xs, ys, style, box) {
    g.save();
    if (box) { g.beginPath(); g.rect(box.x, box.y, box.w, box.h); g.clip(); }
    g.strokeStyle = style.colour; g.lineWidth = style.width || 1.5; g.setLineDash(style.dash || []);
    if (style.alpha !== undefined) g.globalAlpha = style.alpha;
    g.beginPath();
    let pen = false;
    for (let i = 0; i < xs.length; i++) {
      const x = map.X(xs[i]), y = map.Y(ys[i]);
      if (!Number.isFinite(x) || !Number.isFinite(y)) { pen = false; continue; }
      if (pen) g.lineTo(x, y); else { g.moveTo(x, y); pen = true; }
    }
    g.stroke();
    g.restore();
  },

  // Legend entries [{label, colour, dash}] at (x, y).
  legend(g, x, y, entries, u) {
    g.save();
    g.font = `${Math.round(11 * u)}px system-ui, sans-serif`;
    g.textBaseline = 'middle';
    entries.forEach((e, i) => {
      const yy = y + i * 16 * u;
      g.strokeStyle = e.colour; g.fillStyle = e.colour; g.lineWidth = 2 * u; g.setLineDash(e.dash || []);
      if (e.dot) { g.beginPath(); g.arc(x + 9 * u, yy, 3 * u, 0, 2 * Math.PI); g.fill(); }
      else { g.beginPath(); g.moveTo(x, yy); g.lineTo(x + 18 * u, yy); g.stroke(); }
      g.setLineDash([]);
      g.fillStyle = '#c9d1e0';
      g.fillText(e.label, x + 24 * u, yy);
    });
    g.restore();
  },
};
