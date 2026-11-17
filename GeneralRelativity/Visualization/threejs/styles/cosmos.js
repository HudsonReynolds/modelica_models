// The expanding universe in 3D (Examples.CosmicExpansion): a cloud of comoving galaxies whose distances grow
// with the scale factor a(t), with the Hubble sphere, the particle horizon (the observable universe), the
// event horizon and the sphere the cosmic microwave background we see today comes from. Galaxies are
// coloured by whether their light has reached us yet and whether light we send can still reach them. An
// inset plots a(t) and the energy budget (radiation, matter, dark energy). The clock can run in log time so
// the early universe is not over in a blink.
(window.GRStyles = window.GRStyles || []).push({
  id: 'cosmos',
  label: 'Expanding universe',
  order: 8,
  description: 'Comoving galaxies carried apart by the expansion, a(t) from the Friedmann equations, with the Hubble sphere, the particle horizon (edge of the observable universe), the event horizon and the last-scattering sphere. Bright galaxies are ones whose light has reached us; dim ones not yet; blue ones are beyond the event horizon: light we send now never reaches them.',
  options: [
    { id: 'frame', label: 'coordinates', type: 'select', value: 0, choices: [{ value: 0, label: 'proper distance' }, { value: 1, label: 'comoving distance' }] },
    { id: 'logTime', label: 'log time clock', type: 'checkbox', value: 1, help: 'Run the clock in log t, so the radiation and matter eras get screen time' },
    { id: 'galaxies', label: 'galaxies', type: 'range', min: 500, max: 20000, step: 500, value: 6000 },
    { id: 'spheres', label: 'horizon spheres', type: 'checkbox', value: 1 },
    { id: 'inset', label: 'a(t) inset', type: 'checkbox', value: 1 },
  ],
  available(data) {
    return data.records.expansion ? null : 'needs an expansion record (Examples.CosmicExpansion)';
  },
  create(ctx) {
    const { THREE, data, options } = ctx;
    const ex = data.records.expansion[0], cos = (data.records.cosmology || [])[0] || {};
    const n = ex.t.length, T = ex.t, A = ex.a, ETA = ex.eta;
    const interp = (xs, ys, x) => {
      if (x <= xs[0]) return ys[0];
      if (x >= xs[xs.length - 1]) return ys[ys.length - 1];
      let lo = 0, hi = xs.length - 1;
      while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (xs[mid] <= x) lo = mid; else hi = mid; }
      return ys[lo] + (x - xs[lo]) / (xs[hi] - xs[lo]) * (ys[hi] - ys[lo]);
    };
    const logT = T.map(Math.log10), logA = A.map(Math.log10);
    // log-log interpolation of a(t) (a is a power law in t piecewise).
    const aOf = (t) => Math.pow(10, interp(logT, logA, Math.log10(Math.max(t, T[0]))));
    const etaOf = (t) => Math.pow(10, interp(logT, ETA.map((e) => Math.log10(Math.max(e, 1e-30))), Math.log10(Math.max(t, T[0]))));
    const age = cos.age || interp(A, T, 1);
    const eta0 = etaOf(age);
    const mid = Math.floor(n / 2), etaInf = ETA[mid] + ex.dEvent[mid] / A[mid];
    const chiStar = cos.chiStar || (eta0 - interp(A, ETA, 1 / 1091));
    const timeOfClock = (clock) => (options.logTime ? Math.pow(10, clock) : clock);
    const group = new THREE.Group();

    // Galaxies: uniform in a comoving ball a little larger than the event horizon at t -> infinity.
    const R = 1.08 * etaInf, N = Math.round(options.galaxies);
    let seed = 12345;
    const rand = () => { seed = (seed * 1664525 + 1013904223) >>> 0; return seed / 4294967296; };
    const pos = new Float32Array(3 * N), col = new Float32Array(3 * N), chi = new Float32Array(N);
    for (let i = 0; i < N; i++) {
      const r = R * Math.cbrt(rand()), ct = 2 * rand() - 1, ph = 2 * Math.PI * rand(), st = Math.sqrt(1 - ct * ct);
      pos[3 * i] = r * st * Math.cos(ph); pos[3 * i + 1] = r * ct; pos[3 * i + 2] = r * st * Math.sin(ph);
      chi[i] = r;
    }
    const geom = new THREE.BufferGeometry();
    geom.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    geom.setAttribute('color', new THREE.BufferAttribute(col, 3));
    const cloud = new THREE.Points(geom, new THREE.PointsMaterial({ size: 2.2, sizeAttenuation: false, vertexColors: true, transparent: true, opacity: 0.95, depthWrite: false }));
    group.add(cloud);
    group.add(new THREE.Mesh(new THREE.SphereGeometry(0.6, 16, 12), new THREE.MeshBasicMaterial({ color: 0xffffff })));

    // Horizon spheres: a faint shell and a bright equator ring, scaled every frame.
    const sphereSpecs = [
      { key: 'hubble', label: 'Hubble sphere', colour: '#ff6b8a' },
      { key: 'particle', label: 'particle horizon (observable universe)', colour: '#7ee0a1' },
      { key: 'event', label: 'event horizon', colour: '#ffb454' },
      { key: 'lss', label: 'last scattering (the CMB we see today)', colour: '#5fb3ff' },
    ];
    const ringPoints = [];
    for (let k = 0; k <= 128; k++) { const a = 2 * Math.PI * k / 128; ringPoints.push(new THREE.Vector3(Math.cos(a), 0, Math.sin(a))); }
    const spheres = sphereSpecs.map((s) => {
      const g2 = new THREE.Group();
      const colour = new THREE.Color(s.colour);
      g2.add(new THREE.Mesh(new THREE.SphereGeometry(1, 48, 24), new THREE.MeshBasicMaterial({ color: colour, transparent: true, opacity: 0.05, depthWrite: false, side: THREE.DoubleSide })));
      g2.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(ringPoints), new THREE.LineBasicMaterial({ color: colour })));
      const tilted = new THREE.Line(new THREE.BufferGeometry().setFromPoints(ringPoints), new THREE.LineBasicMaterial({ color: colour, transparent: true, opacity: 0.35 }));
      tilted.rotation.x = Math.PI / 2;
      g2.add(tilted);
      group.add(g2);
      return { ...s, obj: g2 };
    });

    const layer = GRPlot.layer(ctx);
    const omegaR = cos.omegaR || 0, omegaM = cos.omegaM || 0.31, omegaL = cos.omegaL || 0.69;
    const fractions = (a) => {
      const r = omegaR / a ** 4, m = omegaM / a ** 3, l = omegaL, s = r + m + l;
      return [r / s, m / s, l / s];
    };
    const epochs = [
      [cos.tEq, 'matter-radiation equality'], [cos.tStar, 'last scattering (CMB)'], [cos.tAcc, 'acceleration begins'], [age, 'today'],
    ].filter(([t]) => t > 0);
    ctx.setNotes([
      `Planck 2018 LCDM: age ${ctx.fmt(age)} Gyr. Today: Hubble sphere ${ctx.fmt(interp(T, ex.dHubble, age))} Gly, observable universe (particle horizon) ${ctx.fmt(eta0)} Gly, event horizon ${ctx.fmt(etaInf - eta0)} Gly, last-scattering sphere ${ctx.fmt(chiStar)} Gly (all proper distances today).`,
      'In comoving coordinates the galaxies stay put and the horizons move: the particle horizon only grows, the comoving Hubble sphere grows while the expansion decelerates and shrinks once it accelerates.',
    ].join('\n'));

    let state = { t: age, a: 1, eta: eta0 };
    function update(clock) {
      const t = timeOfClock(clock), a = aOf(t), eta = etaOf(t);
      state = { t, a, eta };
      const proper = Math.round(options.frame) === 0, s = proper ? a : 1;
      cloud.scale.setScalar(s);
      // Seen (light has reached us), not yet seen, and beyond the event horizon (unreachable from now on).
      const seen = new THREE.Color('#ffd9a0'), unseen = new THREE.Color('#566176'), lost = new THREE.Color('#3b6fd0');
      for (let i = 0; i < N; i++) {
        const c = chi[i] < eta ? seen : chi[i] > etaInf - eta ? lost : unseen;
        col[3 * i] = c.r; col[3 * i + 1] = c.g; col[3 * i + 2] = c.b;
      }
      geom.attributes.color.needsUpdate = true;
      const radii = { hubble: interp(T, ex.dHubble, t) / (proper ? 1 : a), particle: eta, event: etaInf - eta, lss: chiStar };
      for (const sp of spheres) {
        const rr = Math.max(radii[sp.key] * (proper && sp.key !== 'hubble' ? a : 1), 1e-9);
        sp.obj.scale.setScalar(rr);
        sp.obj.visible = !!options.spheres && !(sp.key === 'lss' && t < (cos.tStar || 0));
      }
    }

    function drawInset(r) {
      const { w, h, u } = layer.begin(r);
      const g = layer.g;
      g.clearRect(0, 0, w, h);
      // Legend (top left).
      if (options.spheres) {
        g.fillStyle = 'rgba(11,14,20,0.55)'; g.fillRect(8 * u, 8 * u, 300 * u, 92 * u);
        GRPlot.legend(g, 16 * u, 20 * u, sphereSpecs.map((s) => ({ label: s.label, colour: s.colour })), u);
        g.save(); g.font = `${Math.round(11 * u)}px system-ui`; g.fillStyle = '#c9d1e0';
        g.fillText(`${Math.round(options.frame) === 0 ? 'proper' : 'comoving'} distances; galaxies: seen / not yet / beyond the event horizon`, 16 * u, 92 * u);
        g.restore();
      }
      if (!options.inset) return;
      const W = Math.min(0.42 * w, 520 * u), H = Math.min(0.42 * h, 330 * u), x0 = 12 * u, y0 = h - H - 12 * u;
      g.fillStyle = 'rgba(11,14,20,0.82)'; g.fillRect(x0, y0, W, H);
      const t0 = Math.max(T[0], 1e-7), t1 = T[n - 1];
      const box1 = { x: x0 + 52 * u, y: y0 + 10 * u, w: W - 64 * u, h: 0.52 * H - 26 * u };
      const box2 = { x: box1.x, y: box1.y + box1.h + 12 * u, w: box1.w, h: 0.48 * H - 44 * u };
      const m1 = GRPlot.axes(g, box1, { lo: t0, hi: t1, log: true, noLabels: true, ticks: GRPlot.logTicks(t0, t1).filter((v, i, arr) => i % 2 === 0 || arr.length < 6) },
        { lo: A[0] > 1e-9 ? A[0] : 1e-9, hi: Math.max(A[n - 1], 2), log: true, label: 'a', ticks: [1e-6, 1e-3, 1] }, u);
      GRPlot.line(g, m1, T, A, { colour: '#e6e9ef', width: 1.6 * u }, box1);
      const fr = A.map(fractions);
      const m2 = GRPlot.axes(g, box2, { lo: t0, hi: t1, log: true, label: 'cosmic time (Gyr)', ticks: GRPlot.logTicks(t0, t1).filter((v, i, arr) => i % 2 === 0 || arr.length < 6) },
        { lo: 0, hi: 1, ticks: [0, 1], label: 'share' }, u);
      [['radiation', '#ff8a5c', 0], ['matter', '#5fb3ff', 1], ['dark energy', '#c792ea', 2]].forEach(([label, colour, k]) => {
        GRPlot.line(g, m2, T, fr.map((f) => f[k]), { colour, width: 1.6 * u }, box2);
      });
      GRPlot.legend(g, box2.x + box2.w - 110 * u, box2.y + 10 * u, [{ label: 'radiation', colour: '#ff8a5c' }, { label: 'matter', colour: '#5fb3ff' }, { label: 'dark energy', colour: '#c792ea' }], u * 0.9);
      // Epochs and the clock.
      g.save();
      g.font = `${Math.round(9.5 * u)}px system-ui`;
      for (const [t, label] of epochs) {
        const x = m1.X(t);
        g.strokeStyle = 'rgba(255,255,255,0.25)'; g.setLineDash([2 * u, 3 * u]);
        g.beginPath(); g.moveTo(x, box1.y); g.lineTo(x, box2.y + box2.h); g.stroke();
        g.save(); g.translate(x + 3 * u, box1.y + 4 * u); g.rotate(Math.PI / 2); g.fillStyle = '#8a93a6'; g.fillText(label, 0, 0); g.restore();
      }
      g.setLineDash([]);
      const x = m1.X(Math.max(state.t, t0));
      g.strokeStyle = '#5fb3ff'; g.lineWidth = 2 * u;
      g.beginPath(); g.moveTo(x, box1.y); g.lineTo(x, box2.y + box2.h); g.stroke();
      g.restore();
    }

    return {
      group,
      hidesBodies: false,
      bodyRadius: () => 0.9,
      timeOf: (sample) => (options.logTime ? Math.log10(Math.max(sample.x[0] / (data.header.c || 1), T[0])) : sample.x[0] / (data.header.c || 1)),
      get clockLabel() { return options.logTime ? 'log t' : 't'; },
      place(x) {
        const p = ctx.planeXYZ(x);
        const s = Math.round(options.frame) === 0 ? 1 : 1 / aOf(Math.max(x[0], T[0]));
        return ctx.toWorld(p.X * s, p.Y * s, p.Z * s);
      },
      update(clock) { update(clock); },
      render(r, camera, scene3) {
        r.render(scene3, camera);
        drawInset(r);
        layer.render(r, true);
      },
      camera() {
        const d = 1.9 * eta0;
        return { position: new THREE.Vector3(0.55 * d, 0.45 * d, 0.9 * d), target: new THREE.Vector3(0, 0, 0) };
      },
      readout() {
        const { t, a, eta } = state;
        const [fr, fm, fl] = fractions(a);
        const dom = fr > fm && fr > fl ? 'radiation' : fm > fl ? 'matter' : 'dark energy';
        const H = interp(T, ex.H, t), q = interp(T, ex.q || ex.t.map(() => 0), t);
        return [
          `t = ${t < 1e-3 ? ctx.fmt(t * 1e9) + ' yr' : t < 1 ? ctx.fmt(t * 1e3) + ' Myr' : ctx.fmt(t) + ' Gyr'}, a = ${ctx.fmt(a)}, z = ${ctx.fmt(1 / a - 1)}`,
          `T = ${ctx.fmt((cos.Tcmb || 2.7255) / a)} K, H = ${ctx.fmt(H)} km/s/Mpc, q = ${ctx.fmt(q)} (${q < 0 ? 'accelerating' : 'decelerating'}), ${dom} dominated`,
          `proper distances: Hubble sphere ${ctx.fmt(interp(T, ex.dHubble, t))} Gly, particle horizon ${ctx.fmt(a * eta)} Gly, event horizon ${ctx.fmt(a * (etaInf - eta))} Gly`,
        ];
      },
      onOption(id) {
        if (id === 'logTime' || id === 'galaxies' || id === 'frame') return false;   // rebuild (time mapping, cloud, track placement)
        return true;
      },
      dispose() { layer.dispose(); },
    };
  },
});
