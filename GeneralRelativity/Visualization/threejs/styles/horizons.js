// Light cones and horizons of an expanding universe (after Davis & Lineweaver, PASA 21, 97 (2004)): our
// past light cone, the Hubble sphere, the particle and event horizons and comoving galaxies, in proper
// distance against cosmic time, comoving distance against cosmic time, or comoving distance against
// conformal time (where light moves at 45 degrees and the whole infinite future fits in a finite height).
// Needs the `expansion` record of Examples.CosmicExpansion (lengths in Gly, times in Gyr, c = 1).
(window.GRStyles = window.GRStyles || []).push({
  id: 'horizons',
  label: 'Light cones and horizons',
  order: 9,
  description: 'Spacetime diagram of the expanding universe: our past light cone (the teardrop), the Hubble sphere, the particle horizon (what we can see) and the event horizon (what we can ever reach), with comoving galaxies; in proper or comoving distance, against cosmic or conformal time. The highlighted cone is the light cone of an observer at the clock time.',
  options: [
    { id: 'panel', label: 'diagram', type: 'select', value: 0, choices: [
      { value: 0, label: 'proper distance vs time' }, { value: 1, label: 'comoving distance vs time' }, { value: 2, label: 'comoving distance vs conformal time' }] },
    { id: 'tMax', label: 'time axis (Gyr)', type: 'range', min: 5, max: 60, step: 1, value: 30 },
    { id: 'dMax', label: 'distance axis (Gly)', type: 'range', min: 5, max: 120, step: 1, value: 60 },
    { id: 'galaxies', label: 'comoving galaxies', type: 'checkbox', value: 1 },
  ],
  available(data) {
    return data.records.expansion ? null : 'needs an expansion record (Examples.CosmicExpansion)';
  },
  create(ctx) {
    const { THREE, data, options } = ctx;
    const ex = data.records.expansion[0], cos = (data.records.cosmology || [])[0] || {};
    const n = ex.t.length;
    const T = ex.t, A = ex.a, ETA = ex.eta;
    // Conformal time today and at infinity (from the event horizon a (eta_inf - eta)).
    const interp = (xs, ys, x) => {
      if (x <= xs[0]) return ys[0];
      if (x >= xs[xs.length - 1]) return ys[ys.length - 1];
      let lo = 0, hi = xs.length - 1;
      while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (xs[mid] <= x) lo = mid; else hi = mid; }
      const f = (x - xs[lo]) / (xs[hi] - xs[lo]);
      return ys[lo] + f * (ys[hi] - ys[lo]);
    };
    const age = cos.age || interp(A, T, 1);
    const eta0 = interp(T, ETA, age);
    const etaInf = ETA[Math.floor(n / 2)] + ex.dEvent[Math.floor(n / 2)] / A[Math.floor(n / 2)];
    const etaOfA = (a) => interp(A, ETA, a);
    // Comoving distances of the labelled galaxies: their light reaches us today from redshift z.
    const sources = [[0.5, '#ffd27f'], [1, '#ffb454'], [3, '#ff8a5c'], [10, '#e8697a'], [cos.zStar || 1090, '#5fb3ff']]
      .map(([z, colour]) => ({ z, colour, chi: eta0 - etaOfA(1 / (1 + z)) }));
    const layer = GRPlot.layer(ctx);
    const group = new THREE.Group();
    ctx.setNotes([
      `Planck 2018 LCDM: age ${ctx.fmt(age)} Gyr; today the Hubble sphere is ${ctx.fmt(cos.hubbleRadiusToday || interp(T, ex.dHubble, age))} Gly, the particle horizon ${ctx.fmt(eta0)} Gly and the event horizon ${ctx.fmt(etaInf - eta0)} Gly away.`,
      'Galaxies beyond the Hubble sphere recede faster than light, yet we see some of them: light crosses into the sphere as it grows. Light emitted today beyond the event horizon will never reach us.',
    ].join('\n'));

    let lastKey = '';
    function draw(r) {
      const { w, h, u } = layer.begin(r);
      const t = ctx.clock;
      const key = [w, h, t, JSON.stringify(options)].join('|');
      if (key === lastKey) return;
      lastKey = key;
      const g = layer.g;
      g.fillStyle = '#0b0e14';
      g.fillRect(0, 0, w, h);
      const panel = Math.round(options.panel);
      const conformal = panel === 2, proper = panel === 0;
      const box = { x: 70 * u, y: 34 * u, w: w - 150 * u, h: h - 90 * u };
      const yMax = conformal ? Math.min(etaInf * 1.02, options.tMax * 2.2) : options.tMax;
      const map = GRPlot.axes(g, box,
        { lo: -options.dMax, hi: options.dMax, label: proper ? 'proper distance (Gly)' : 'comoving distance (Gly)' },
        { lo: 0, hi: yMax, label: conformal ? 'conformal time (Gyr)' : 'cosmic time (Gyr)' }, u);
      // Rows of the table up to the top of the plot, as (time-axis value, a, eta).
      const rows = [];
      for (let i = 0; i < n; i++) {
        const y = conformal ? ETA[i] : T[i];
        rows.push({ y, a: A[i], eta: ETA[i], t: T[i], hub: ex.dHubble[i] });
        if (y > yMax) break;
      }
      const ys = rows.map((p) => p.y);
      const dist = (p, chi) => (proper ? p.a * chi : chi);   // a comoving distance shown on the chosen axis
      const curve = (fn, style, mirror = true, filter = () => true) => {
        const sel = rows.filter(filter);
        const xs = sel.map(fn), yv = sel.map((p) => p.y);
        GRPlot.line(g, map, xs, yv, style, box);
        if (mirror) GRPlot.line(g, map, xs.map((x) => -x), yv, style, box);
      };
      // Comoving galaxies.
      if (options.galaxies) {
        for (const s of sources) curve((p) => dist(p, s.chi), { colour: s.colour, width: 1 * u, alpha: 0.7 });
        for (let chi = 10; chi <= 3 * options.dMax; chi += 10) curve((p) => dist(p, chi), { colour: '#3a4458', width: 0.8 * u });
      }
      // Horizons.
      curve((p) => (proper ? p.hub : p.hub / p.a), { colour: '#ff6b8a', width: 2 * u });
      curve((p) => dist(p, p.eta), { colour: '#7ee0a1', width: 2 * u, dash: [8 * u, 5 * u] });
      curve((p) => dist(p, etaInf - p.eta), { colour: '#ffb454', width: 2 * u, dash: [3 * u, 4 * u] });
      // Our light cone today (past teardrop and future cone).
      curve((p) => dist(p, Math.abs(eta0 - p.eta)), { colour: '#e6e9ef', width: 2 * u });
      // The light cone of an observer at the clock time.
      const etaNow = interp(T, ETA, t);
      curve((p) => dist(p, etaNow - p.eta), { colour: '#5fb3ff', width: 2.5 * u }, true, (p) => p.eta <= etaNow);
      // Now and today lines, last scattering.
      const yNow = conformal ? etaNow : t, yToday = conformal ? eta0 : age;
      g.save();
      g.strokeStyle = '#5fb3ff'; g.lineWidth = 1.5 * u; g.setLineDash([6 * u, 4 * u]);
      g.beginPath(); g.moveTo(box.x, map.Y(yNow)); g.lineTo(box.x + box.w, map.Y(yNow)); g.stroke();
      g.strokeStyle = 'rgba(255,255,255,0.4)'; g.setLineDash([2 * u, 4 * u]);
      g.beginPath(); g.moveTo(box.x, map.Y(yToday)); g.lineTo(box.x + box.w, map.Y(yToday)); g.stroke();
      g.restore();
      // Right axis: the scale factor.
      g.save();
      g.font = `${Math.round(10 * u)}px system-ui`; g.fillStyle = '#8a93a6'; g.textBaseline = 'middle';
      for (const a of [0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5, 10]) {
        const y = conformal ? etaOfA(a) : interp(A, T, a);
        if (y > yMax || y < 0) continue;
        g.fillText(`a = ${a}`, box.x + box.w + 6 * u, map.Y(y));
      }
      g.fillText('scale factor', box.x + box.w + 6 * u, box.y - 10 * u);
      g.restore();
      // Labels for the galaxies along the top of their lines.
      if (options.galaxies) {
        g.save();
        g.font = `${Math.round(10 * u)}px system-ui`; g.textAlign = 'center';
        for (const s of sources) {
          const p = rows[rows.length - 1];
          const x = map.X(dist(p, s.chi));
          if (x < box.x + box.w - 10 * u && x > box.x) { g.fillStyle = s.colour; g.fillText(`z = ${s.z >= 1000 ? Math.round(s.z) : s.z}`, x, box.y - 6 * u); }
        }
        g.restore();
      }
      GRPlot.legend(g, box.x + 10 * u, box.y + 14 * u, [
        { label: 'our light cone today', colour: '#e6e9ef' },
        { label: `light cone at t = ${ctx.fmt(t)} Gyr`, colour: '#5fb3ff' },
        { label: 'Hubble sphere (recession at c)', colour: '#ff6b8a' },
        { label: 'particle horizon', colour: '#7ee0a1', dash: [8 * u, 5 * u] },
        { label: 'event horizon', colour: '#ffb454', dash: [3 * u, 4 * u] },
        { label: 'comoving galaxies (z of the light we get today)', colour: '#ffb454' },
      ], u);
      g.save();
      g.font = `${Math.round(13 * u)}px system-ui`; g.fillStyle = '#e6e9ef'; g.textBaseline = 'top';
      g.fillText(['Proper distance', 'Comoving distance', 'Comoving distance, conformal time'][panel], box.x, 8 * u);
      g.restore();
    }
    return {
      group,
      hidesTracks: true,
      hidesBodies: true,
      place(x) { return ctx.toWorld(x[1], x[2], x[3]); },
      render(r) { draw(r); layer.render(r); },
      readout() {
        const t = ctx.clock, a = interp(T, A, t), eta = interp(T, ETA, t);
        return [`t = ${ctx.fmt(t)} Gyr, a = ${ctx.fmt(a)}, z = ${ctx.fmt(1 / a - 1)}`,
          `Hubble sphere ${ctx.fmt(interp(T, ex.dHubble, t))} Gly, particle horizon ${ctx.fmt(a * eta)} Gly, event horizon ${ctx.fmt(a * (etaInf - eta))} Gly (proper)`];
      },
      onOption() { lastKey = ''; return true; },
      dispose() { layer.dispose(); },
    };
  },
});
