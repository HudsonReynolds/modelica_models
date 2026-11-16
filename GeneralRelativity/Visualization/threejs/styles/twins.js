// Twin paradox style: a Minkowski spacetime diagram (x across, ct up, light at 45 degrees) of a twin who
// stays home and a twin who flies out and back with constant proper acceleration a, cruise speed beta and
// turnaround distance D. Everything is recomputed live from the sliders with the closed-form hyperbolic
// motion that Examples.TwinParadox integrates and asserts in Modelica: per phase of proper acceleration
// alpha, t and x advance by (sinh, cosh)(eta0 + alpha s) differences over alpha, eta = rapidity.
// Units: years and light-years, c = 1 (1 g = 1.0323 ly/yr^2).
(function () {
  const G_LY_PER_YR2 = 9.80665 * Math.pow(365.25 * 86400, 2) / 9.4607304725808e15;

  function advance(t, x, eta, alpha, s) {
    if (alpha !== 0) {
      return [t + (Math.sinh(eta + alpha * s) - Math.sinh(eta)) / alpha, x + (Math.cosh(eta + alpha * s) - Math.cosh(eta)) / alpha, eta + alpha * s];
    }
    return [t + Math.cosh(eta) * s, x + Math.sinh(eta) * s, eta];
  }
  // Out at +a to the cruise rapidity, cruise, -a for twice as long (stop and come back), cruise, +a to rest.
  function makeTrip(beta, D, a) {
    const etaCruise = Math.atanh(Math.min(beta, 0.9999));
    const etaMax = Math.min(etaCruise, Math.acosh(1 + a * D / 2));
    const tauA = etaMax / a;
    const xa = (Math.cosh(etaMax) - 1) / a;
    const bc = Math.tanh(etaMax);
    const tc = Math.max(0, (D - 2 * xa) / Math.max(bc, 1e-12));
    const tauC = tc / Math.cosh(etaMax);
    const spec = [[a, tauA], [0, tauC], [-a, 2 * tauA], [0, tauC], [a, tauA]];
    const phases = [];
    let t = 0, x = 0, eta = 0, tau = 0;
    for (const [alpha, dtau] of spec) {
      phases.push({ t, x, eta, tau, alpha, dtau });
      [t, x, eta] = advance(t, x, eta, alpha, dtau);
      tau += dtau;
    }
    return { phases, T: t, tauTotal: tau, etaMax, betaReached: bc, reached: etaMax >= etaCruise - 1e-12, xEnd: x, a, D, beta };
  }
  function atTau(trip, tau) {
    tau = Math.min(Math.max(tau, 0), trip.tauTotal);
    let p = trip.phases[0];
    for (const q of trip.phases) if (q.tau <= tau) p = q;
    const [t, x, eta] = advance(p.t, p.x, p.eta, p.alpha, tau - p.tau);
    return { t, x, eta, tau };
  }
  function atT(trip, t) {
    // t(tau) is increasing: bisection.
    let lo = 0, hi = trip.tauTotal;
    if (t >= trip.T) return atTau(trip, hi);
    for (let k = 0; k < 64; k++) { const mid = 0.5 * (lo + hi); if (atTau(trip, mid).t < t) lo = mid; else hi = mid; }
    return atTau(trip, 0.5 * (lo + hi));
  }
  // Light from home emitted at t_e (moving +x) meets the traveller when t - t_e = x(t).
  function homeSignalArrival(trip, te) {
    let lo = te, hi = trip.T;
    const f = (t) => (t - te) - atT(trip, t).x;
    if (f(hi) < 0) return null;
    for (let k = 0; k < 50; k++) { const mid = 0.5 * (lo + hi); if (f(mid) < 0) lo = mid; else hi = mid; }
    return 0.5 * (lo + hi);
  }

  window.GRStyles = window.GRStyles || [];
  window.GRStyles.push({
    id: 'twins',
    label: 'Twin paradox (Minkowski)',
    order: 6,
    description: 'Minkowski diagram of the twin paradox: x across, ct up, light at 45 degrees. Move the sliders: the traveller\'s worldline, proper-time ticks, yearly light signals and lines of simultaneity are recomputed at once.',
    options: [
      { id: 'beta', label: 'cruise speed v/c', type: 'range', min: 0.05, max: 0.995, step: 0.005, value: 0.8 },
      { id: 'distance', label: 'turnaround (ly)', type: 'range', min: 0.2, max: 20, step: 0.1, value: 4 },
      { id: 'acceleration', label: 'acceleration (g)', type: 'range', min: 0.1, max: 10, step: 0.05, value: 1, help: 'Proper acceleration in units of g = 9.80665 m/s^2 = 1.032 ly/yr^2' },
      { id: 'tick', label: 'tick every (yr)', type: 'range', min: 0.25, max: 5, step: 0.25, value: 1 },
      { id: 'signals', label: 'birthday signals', type: 'checkbox', value: 1, help: 'Each twin sends a light flash at every tick of their own clock' },
      { id: 'simultaneity', label: 'simultaneity', type: 'checkbox', value: 1, help: 'The traveller\'s line of simultaneity now (slope v/c); it swings during the turnaround' },
      { id: 'simulated', label: 'Modelica tracks', type: 'checkbox', value: 1, help: 'Overlay the worldlines integrated by Examples.TwinParadox (scene parameters)' },
    ],
    available(data) {
      return data.view && data.view.style === 'twins' ? null : 'needs a twin-paradox scene (Examples.TwinParadox)';
    },
    create(ctx) {
      const { THREE, data, options } = ctx;
      const group = new THREE.Group();
      const specs = data.trackSpecs;
      const homeSpec = specs.find(s => /home/i.test(s.name)) || specs[0];
      const travSpec = specs.find(s => s !== homeSpec) || specs[0];
      const homeColour = new THREE.Color(homeSpec && homeSpec.color || '#5fb3ff');
      const travColour = new THREE.Color(travSpec && travSpec.color || '#ffb454');
      const sceneOptions = (data.view && data.view.options) || {};
      let trip = null;
      const bodyR = () => 0.012 * Math.max(trip.T, trip.D + 2);
      const homeBody = new THREE.Mesh(new THREE.SphereGeometry(1, 24, 16), new THREE.MeshBasicMaterial({ color: homeColour }));
      const travBody = new THREE.Mesh(new THREE.SphereGeometry(1, 24, 16), new THREE.MeshBasicMaterial({ color: travColour }));
      const nowLine = new THREE.Line(new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(), new THREE.Vector3()]),
        new THREE.LineDashedMaterial({ color: 0xffffff, dashSize: 0.15, gapSize: 0.1, transparent: true, opacity: 0.5 }));
      const simLine = new THREE.Line(new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(), new THREE.Vector3()]),
        new THREE.LineBasicMaterial({ color: travColour, transparent: true, opacity: 0.85 }));
      const coneLines = new THREE.LineSegments(new THREE.BufferGeometry().setFromPoints([0, 1, 2, 3].map(() => new THREE.Vector3())),
        new THREE.LineBasicMaterial({ color: 0xf9e2af, transparent: true, opacity: 0.8 }));
      let staticGroup = new THREE.Group();
      group.add(staticGroup, homeBody, travBody, nowLine, simLine, coneLines);

      const V = (x, t) => new THREE.Vector3(x, t, 0);
      // Text labels as sprites drawn on a canvas.
      const label = (text, colour, x, t, size) => {
        const canvas = document.createElement('canvas');
        canvas.width = 512; canvas.height = 96;
        const c2 = canvas.getContext('2d');
        c2.font = '600 52px system-ui, sans-serif';
        c2.fillStyle = colour;
        c2.textBaseline = 'middle';
        c2.fillText(text, 8, 48);
        const sprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: new THREE.CanvasTexture(canvas), transparent: true, depthWrite: false }));
        sprite.scale.set(size * 512 / 96, size, 1);
        sprite.center.set(0, 0.5);
        sprite.position.set(x, t, 0);
        return sprite;
      };
      const segs = (points, colour, opacity = 1) => {
        const g = new THREE.BufferGeometry().setFromPoints(points);
        return new THREE.LineSegments(g, new THREE.LineBasicMaterial({ color: colour, transparent: opacity < 1, opacity }));
      };
      let signalStats = null;
      let built = false;
      const framing = () => {
        const W = trip.D + 2, H = trip.T;
        const d = 0.62 * Math.max(H, W * 0.75) / Math.tan(THREE.MathUtils.degToRad(25));
        return { position: new THREE.Vector3(trip.D / 2, H / 2, d), target: new THREE.Vector3(trip.D / 2, H / 2, 0) };
      };

      function build() {
        group.remove(staticGroup);
        staticGroup.traverse(o => { if (o.geometry) o.geometry.dispose(); if (o.material) { if (o.material.map) o.material.map.dispose(); o.material.dispose(); } });
        staticGroup = new THREE.Group();
        group.add(staticGroup);
        trip = makeTrip(options.beta, options.distance, options.acceleration * G_LY_PER_YR2);
        const T = trip.T, D = trip.D;
        const xLo = -1, xHi = D + 1;
        // Grid: one line per year / light-year (coarser for long trips).
        const step = Math.max(T, D) > 40 ? 5 : 1;
        const grid = [];
        for (let x = Math.ceil(xLo / step) * step; x <= xHi; x += step) grid.push(V(x, 0), V(x, T));
        for (let t = 0; t <= T; t += step) grid.push(V(xLo, t), V(xHi, t));
        staticGroup.add(segs(grid, 0x2a3346));
        staticGroup.add(segs([V(xLo, 0), V(xHi, 0), V(0, 0), V(0, T * 1.04)], 0x8a93a6));
        // Light cone of the departure event.
        staticGroup.add(segs([V(0, 0), V(Math.min(xHi, T), Math.min(xHi, T)), V(0, 0), V(Math.max(xLo, -T), -Math.max(xLo, -T))], 0xf9e2af, 0.35));
        // Worldlines.
        staticGroup.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints([V(0, 0), V(0, T)]), new THREE.LineBasicMaterial({ color: homeColour })));
        const path = [];
        for (let k = 0; k <= 600; k++) { const e = atTau(trip, trip.tauTotal * k / 600); path.push(V(e.x, e.t)); }
        staticGroup.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(path), new THREE.LineBasicMaterial({ color: travColour })));
        // Proper-time ticks: a short bar across each worldline at every tick of its own clock.
        const w = 0.012 * Math.max(T, D + 2);
        const tickPts = [], tick = options.tick;
        const nHome = Math.floor(T / tick + 1e-9), nTrav = Math.floor(trip.tauTotal / tick + 1e-9);
        for (let k = 1; k <= nHome; k++) tickPts.push(V(-w, k * tick), V(w, k * tick));
        staticGroup.add(segs(tickPts, homeColour));
        const travTicks = [];
        for (let k = 1; k <= nTrav; k++) { const e = atTau(trip, k * tick); travTicks.push(V(e.x - w, e.t), V(e.x + w, e.t)); }
        staticGroup.add(segs(travTicks, travColour));
        // Birthday signals at 45 degrees.
        signalStats = { home: [], trav: [] };
        if (options.signals && nHome + nTrav <= 400) {
          const fromHome = [], fromTrav = [];
          for (let k = 1; k <= nHome; k++) {
            const te = k * tick, tr = homeSignalArrival(trip, te);
            if (tr === null) continue;
            const e = atT(trip, tr);
            fromHome.push(V(0, te), V(e.x, tr));
            signalStats.home.push({ sent: te, received: tr, travellerAge: e.tau });
          }
          for (let k = 1; k <= nTrav; k++) {
            const e = atTau(trip, k * tick);
            fromTrav.push(V(e.x, e.t), V(0, e.t + e.x));
            signalStats.trav.push({ sent: k * tick, received: e.t + e.x });
          }
          staticGroup.add(segs(fromHome, homeColour, 0.45), segs(fromTrav, travColour, 0.45));
        }
        const ls = 0.028 * Math.max(T, D + 2);
        staticGroup.add(label('x (light-years)', '#8a93a6', xHi - 4.2 * ls, -1.2 * ls, ls),
          label('ct (years)', '#8a93a6', 0.4 * ls, T * 1.04, ls),
          label('home', '#' + homeColour.getHexString(), -0.2 * ls - 3.0 * ls, 0.5 * T, ls),
          label('traveller', '#' + travColour.getHexString(), D + 0.3 * ls, atTau(trip, 0.5 * trip.tauTotal).t + 1.5 * ls, ls),
          label('light', '#f9e2af', Math.min(xHi, T) * 0.55 + 0.4 * ls, Math.min(xHi, T) * 0.55, 0.8 * ls));
        ctx.setTimeRange(0, T);
        const notes = [];
        if (!trip.reached) notes.push(`Turnaround too close to reach v = ${options.beta.toFixed(3)} c at this acceleration: peak speed ${trip.betaReached.toFixed(4)} c, no cruise.`);
        const s = sceneOptions;
        if (s.beta !== undefined) {
          const same = Math.abs(s.beta - options.beta) < 1e-9 && Math.abs(s.distance - options.distance) < 1e-9 && Math.abs(s.acceleration - options.acceleration) < 1e-9;
          notes.push(same ? 'Sliders match the simulated scene: the Modelica worldlines lie on the live curves.'
            : `The Modelica tracks were simulated for v = ${s.beta} c, D = ${s.distance} ly, a = ${s.acceleration} g.`);
        }
        ctx.setNotes(notes.join('\n'));
        // The diagram is flat: refit the view whenever the sliders change its size.
        if (built) {
          const f = framing();
          ctx.camera.position.copy(f.position);
          ctx.controls.target.copy(f.target);
          ctx.controls.update();
        }
        built = true;
        homeBody.scale.setScalar(bodyR());
        travBody.scale.setScalar(bodyR());
        place(ctx.clock);
      }

      function place(t) {
        const e = atT(trip, t);
        homeBody.position.copy(V(0, t));
        travBody.position.copy(V(e.x, e.t));
        const xLo = -1, xHi = trip.D + 1;
        const nAttr = nowLine.geometry.getAttribute('position');
        nAttr.setXYZ(0, xLo, t, 0); nAttr.setXYZ(1, xHi, t, 0); nAttr.needsUpdate = true;
        nowLine.computeLineDistances();
        // Traveller's simultaneity: t' = const, i.e. t - t_e = (v/c)(x - x_e).
        const b = Math.tanh(e.eta);
        const sAttr = simLine.geometry.getAttribute('position');
        sAttr.setXYZ(0, xLo, e.t + b * (xLo - e.x), 0); sAttr.setXYZ(1, xHi, e.t + b * (xHi - e.x), 0); sAttr.needsUpdate = true;
        simLine.visible = !!options.simultaneity;
        const h = 0.06 * Math.max(trip.T, trip.D + 2);
        const cAttr = coneLines.geometry.getAttribute('position');
        cAttr.setXYZ(0, e.x, e.t, 0); cAttr.setXYZ(1, e.x + h, e.t + h, 0); cAttr.setXYZ(2, e.x, e.t, 0); cAttr.setXYZ(3, e.x - h, e.t + h, 0);
        cAttr.needsUpdate = true;
      }

      build();
      return {
        group,
        hidesBodies: true,
        get hidesTracks() { return !options.simulated; },
        // Modelica samples are {ct, x, y, z}: the same diagram coordinates.
        place(x) { return V(x[1], x[0]); },
        update(clock) { place(clock); },
        camera() { return framing(); },
        readout() {
          const t = ctx.clock, e = atT(trip, t);
          const lines = [
            `peak speed ${trip.betaReached.toFixed(4)} c (gamma ${Math.cosh(trip.etaMax).toFixed(4)}), turnaround ${trip.D.toFixed(2)} ly`,
            `now: home ${t.toFixed(3)} yr, traveller ${e.tau.toFixed(3)} yr, at x = ${e.x.toFixed(3)} ly`,
            `reunion: home ${trip.T.toFixed(4)} yr, traveller ${trip.tauTotal.toFixed(4)} yr`,
            `         the traveller is ${(trip.T - trip.tauTotal).toFixed(4)} yr younger`,
          ];
          if (signalStats && (signalStats.home.length || signalStats.trav.length)) {
            const homeGot = signalStats.trav.filter(sg => sg.received <= t).length;
            const travGot = signalStats.home.filter(sg => sg.received <= t).length;
            lines.push(`signals received so far: home ${homeGot}, traveller ${travGot}`);
          }
          // The Modelica run, read from the scene's traveller track.
          if (travSpec && travSpec.samples.length) {
            const last = travSpec.samples[travSpec.samples.length - 1];
            lines.push(`Modelica run: reunion at t = ${(last.x[0] / (data.header.c || 1)).toFixed(4)} yr, traveller ${last.lambda.toFixed(4)} yr`);
          }
          return lines;
        },
        onOption(id) {
          if (id === 'simulated') { ctx.refreshVisibility(); return true; }
          build();
          return true;
        },
      };
    },
  });
})();
