// Spacetime style: a 2+1 diagram. The equatorial plane is horizontal and ct runs upward on the same scale
// as length (one unit up = one unit of ct), so light far from any mass moves at exactly 45 degrees and a
// worldline can never look faster than its light cone. The view is a window of ct around "now" that
// scrolls as time runs. Cone openings come from the metric (coordinate speeds of light): in
// Schwarzschild coordinates they close up toward the horizon, in Kerr they tilt with frame dragging.
(window.GRStyles = window.GRStyles || []).push({
  id: 'spacetime',
  label: 'Spacetime diagram (2+1)',
  order: 4,
  description: 'Worldlines in a 2+1 spacetime diagram: the equatorial plane horizontal, ct up on the same scale as length (light at 45 degrees). The view scrolls with time. Light cones narrow and tilt where the metric slows and drags light.',
  options: [
    { id: 'window', label: 'time window', type: 'range', min: 0.3, max: 6, step: 0.1, value: 1.5, help: 'Half-height of the visible ct window around now, relative to the view size (the scale stays 1 ct = 1 length)' },
    { id: 'cones', label: 'light cones', type: 'checkbox', value: 1 },
    { id: 'coneSize', label: 'cone size', type: 'range', min: 0.01, max: 0.2, step: 0.005, value: 0.05 },
    { id: 'rings', label: 'cone rings', type: 'range', min: 1, max: 10, step: 1, value: 6 },
    { id: 'bodyCones', label: 'cones on bodies', type: 'checkbox', value: 1, help: 'The future light cone at each body\'s current event' },
    { id: 'now', label: 'now plane', type: 'checkbox', value: 1 },
  ],
  available(data) {
    return data.trackSpecs.length ? null : 'needs at least one trajectory';
  },
  create(ctx) {
    const { THREE, data, options } = ctx;
    const extent = ctx.extent;
    const group = new THREE.Group();
    const W = options.window * extent;               // visible ct above and below now
    const cLight = data.header.c || 1;
    const nowCt = () => cLight * ctx.clock;
    const notes = ['One unit up is one unit of ct, so light far from the mass moves at 45 degrees; the view shows ' + ctx.fmt(W) + ' of ct either side of now and scrolls.'];
    // Clip everything (worldlines included) to the window.
    const renderer = ctx.renderer;
    const previousClipping = renderer.clippingPlanes;
    renderer.clippingPlanes = [new THREE.Plane(new THREE.Vector3(0, -1, 0), W), new THREE.Plane(new THREE.Vector3(0, 1, 0), W)];

    // Time axis with ticks that scroll with ct.
    const tickStep = (() => { const raw = W / 4, p = Math.pow(10, Math.floor(Math.log10(raw))); return [1, 2, 5, 10].map(m => m * p).find(v => v >= raw); })();
    const ticks = new THREE.Group();
    {
      const pts = [new THREE.Vector3(0, -W - tickStep, 0), new THREE.Vector3(0, W + tickStep, 0)];
      for (let k = -Math.ceil(W / tickStep) - 1; k <= Math.ceil(W / tickStep) + 1; k++) {
        pts.push(new THREE.Vector3(-0.02 * extent, k * tickStep, 0), new THREE.Vector3(0.02 * extent, k * tickStep, 0));
      }
      ticks.add(new THREE.LineSegments(new THREE.BufferGeometry().setFromPoints(pts), new THREE.LineBasicMaterial({ color: 0x8a93a6 })));
    }
    group.add(ticks);
    notes.push(`Axis ticks every ${ctx.fmt(tickStep)} of ct.`);
    // Spheres become cylinders (their world tubes) spanning the window.
    for (const s of data.spheres) {
      const color = new THREE.Color(s.color || '#888888');
      const horizon = /horizon/i.test(s.name);
      const cyl = new THREE.Mesh(new THREE.CylinderGeometry(s.radius, s.radius, 2 * W, 64, 1, true),
        new THREE.MeshBasicMaterial({ color: horizon ? 0x000000 : color, transparent: true, opacity: horizon ? 0.75 : 0.18, side: THREE.DoubleSide, depthWrite: false }));
      group.add(cyl);
      const pts = [];
      for (let m = 0; m < 96; m++) {
        const a0 = 2 * Math.PI * m / 96, a1 = 2 * Math.PI * (m + 1) / 96;
        for (const y of [-W, W]) pts.push(new THREE.Vector3(s.radius * Math.cos(a0), y, s.radius * Math.sin(a0)), new THREE.Vector3(s.radius * Math.cos(a1), y, s.radius * Math.sin(a1)));
      }
      for (let m = 0; m < 8; m++) {
        const a = 2 * Math.PI * m / 8;
        pts.push(new THREE.Vector3(s.radius * Math.cos(a), -W, s.radius * Math.sin(a)), new THREE.Vector3(s.radius * Math.cos(a), W, s.radius * Math.sin(a)));
      }
      group.add(new THREE.LineSegments(new THREE.BufferGeometry().setFromPoints(pts),
        new THREE.LineBasicMaterial({ color: horizon ? 0x777777 : color, transparent: true, opacity: 0.45 })));
    }

    // "Now" plane.
    const nowGroup = new THREE.Group();
    {
      const disk = new THREE.Mesh(new THREE.CircleGeometry(extent, 96), new THREE.MeshBasicMaterial({ color: 0x5fb3ff, side: THREE.DoubleSide, transparent: true, opacity: 0.08, depthWrite: false }));
      disk.rotation.x = -Math.PI / 2;
      const rim = [];
      for (let k = 0; k <= 128; k++) rim.push(new THREE.Vector3(extent * Math.cos(2 * Math.PI * k / 128), 0, extent * Math.sin(2 * Math.PI * k / 128)));
      nowGroup.add(disk, new THREE.Line(new THREE.BufferGeometry().setFromPoints(rim), new THREE.LineBasicMaterial({ color: 0x5fb3ff, transparent: true, opacity: 0.5 })));
    }
    nowGroup.visible = !!options.now;
    group.add(nowGroup);

    // Light cones on the now plane and at each body, from the slice's coordinate speeds of light; a flat
    // (Minkowski) scene without a slice gets exact 45-degree cones.
    const coneGroup = new THREE.Group();
    const slice = data.slice;
    const flat = !(slice && slice.fields.coneOut) && data.metric && data.metric.name === 'minkowski';
    let speeds = null;   // (X, Y) -> {so, si, sp, sr} or null
    let rValid = 0;
    if (slice && slice.fields.coneOut) {
      const get = (name) => ctx.sampler(slice, name, 'coneValid');
      const out = get('coneOut'), inw = get('coneIn'), pro = get('conePro'), retro = get('coneRetro');
      speeds = (X, Y) => {
        const so = out(X, Y), si = inw(X, Y), sp = pro(X, Y), sr = retro(X, Y);
        return [so, si, sp, sr].some(v => v === null) ? null : { so, si, sp, sr };
      };
      const rDead = Math.max(slice.rMin || 0, 1e-6 * extent);
      // Innermost radius with valid cones along +x.
      rValid = rDead;
      for (let R = rDead; R < extent; R += extent / 400) if (out(R, 0) !== null) { rValid = R; break; }
    } else if (flat) {
      speeds = () => ({ so: 1, si: -1, sp: 1, sr: -1 });
      notes.push('Flat spacetime: every light cone is exactly 45 degrees.');
    } else {
      notes.push('No slice in the scene: light cones need Visualization.writeSpacetimeSlice.');
    }
    const h = options.coneSize * extent;
    const segs = 32;
    // Line segments of one cone with its apex at (X, Y) on the now plane (height 0 of the cone group).
    const conePoints = (X, Y, sp, y0, pts, cols) => {
      const R = Math.hypot(X, Y);
      const rx = R > 0 ? X / R : 1, ry = R > 0 ? Y / R : 0, tx = -ry, ty = rx;
      const base = [];
      for (let m = 0; m <= segs; m++) {
        const psi = 2 * Math.PI * m / segs, c = Math.cos(psi), sn = Math.sin(psi);
        const offR = c >= 0 ? sp.so * c : -sp.si * c;
        const offT = sn >= 0 ? sp.sp * sn : -sp.sr * sn;
        base.push([X + h * (offR * rx + offT * tx), Y + h * (offR * ry + offT * ty)]);
      }
      const col = ctx.colormap(sp.so);
      const push = (a, ya, b, yb) => { pts.push(a[0], ya, -a[1], b[0], yb, -b[1]); cols.push(col.r, col.g, col.b, col.r, col.g, col.b); };
      for (let m = 0; m < segs; m++) push(base[m], y0 + h, base[m + 1], y0 + h);
      for (let m = 0; m < segs; m += 4) push([X, Y], y0, base[m], y0 + h);
    };
    if (speeds) {
      const rings = Math.round(options.rings);
      // Flat space has no centre to resolve: start the rings further out so they do not pile up at the origin.
      const rLow = Math.max(rValid * 1.08, (flat ? 0.2 : 0.02) * extent), rHigh = 0.85 * extent;
      const pts = [], cols = [];
      let sMin = Infinity;
      for (let k = 0; k < rings; k++) {
        const R = rings === 1 ? 0.5 * (rLow + rHigh) : rLow * Math.pow(rHigh / rLow, k / (rings - 1));
        for (let m = 0; m < 8; m++) {
          const a = 2 * Math.PI * (m + 0.5 * (k % 2)) / 8;
          const X = R * Math.cos(a), Y = R * Math.sin(a), sp = speeds(X, Y);
          if (!sp) continue;
          conePoints(X, Y, sp, 0, pts, cols);
          sMin = Math.min(sMin, sp.so);
        }
      }
      const geom = new THREE.BufferGeometry();
      geom.setAttribute('position', new THREE.Float32BufferAttribute(pts, 3));
      geom.setAttribute('color', new THREE.Float32BufferAttribute(cols, 3));
      coneGroup.add(new THREE.LineSegments(geom, new THREE.LineBasicMaterial({ vertexColors: true })));
      ctx.drawScale(0, 1, 'outward light speed (cone colour)');
      if (!flat && Number.isFinite(sMin)) notes.push(`Slowest outward light shown: ${ctx.fmt(sMin)} c (coordinate speed).`);
    }
    // Cones riding on the bodies (in world coordinates, rebuilt every frame).
    const bodyConeGeom = new THREE.BufferGeometry();
    const bodyCones = new THREE.LineSegments(bodyConeGeom, new THREE.LineBasicMaterial({ vertexColors: true }));
    bodyCones.visible = !!options.cones && !!options.bodyCones && !!speeds;
    group.add(bodyCones);
    function updateBodyCones() {
      if (!bodyCones.visible) return;
      const pts = [], cols = [];
      for (const b of ctx.bodiesNow()) {
        if (b.hidden || !b.position) continue;
        const p = ctx.planeXYZ(b.sample.x), sp = speeds(p.X, p.Y);
        if (sp) conePoints(p.X, p.Y, sp, b.position.y, pts, cols);
      }
      bodyConeGeom.setAttribute('position', new THREE.Float32BufferAttribute(pts, 3));
      bodyConeGeom.setAttribute('color', new THREE.Float32BufferAttribute(cols, 3));
    }
    coneGroup.visible = !!options.cones;
    group.add(coneGroup);
    ctx.setNotes(notes.join('\n'));

    return {
      group,
      // Height = ct relative to now: worldlines scroll down as time runs.
      dynamicPlacement: true,
      place(x) {
        const p = ctx.planeXYZ(x);
        const v = ctx.toWorld(p.X, p.Y, 0);
        v.y = x[0] - nowCt();
        return v;
      },
      update() {
        ticks.position.y = -(((nowCt() % tickStep) + tickStep) % tickStep);
        updateBodyCones();
      },
      camera() {
        return { position: new THREE.Vector3(0.35 * extent, 0.55 * W, 2.3 * Math.max(extent, 0.8 * W)), target: new THREE.Vector3(0, 0, 0) };
      },
      readout() {
        return [`now: ct = ${ctx.fmt(nowCt())}; window +-${ctx.fmt(W)} of ct`];
      },
      dispose() { renderer.clippingPlanes = previousClipping || []; },
      onOption(id, v) {
        if (id === 'cones' || id === 'bodyCones') {
          coneGroup.visible = !!options.cones;
          bodyCones.visible = !!options.cones && !!options.bodyCones && !!speeds;
          return true;
        }
        if (id === 'now') { nowGroup.visible = !!v; return true; }
        return false;
      },
    };
  },
});
