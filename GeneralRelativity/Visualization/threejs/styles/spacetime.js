// Spacetime style: a 2+1 diagram. The equatorial plane is horizontal and coordinate time t runs upward,
// so bodies trace worldlines. Light cones sit on the "now" plane, drawn so that light far from any mass
// moves at 45 degrees; their opening comes from the metric (coordinate speeds of light), so in
// Schwarzschild coordinates they close up toward the horizon and tilt with frame dragging in Kerr.
(window.GRStyles = window.GRStyles || []).push({
  id: 'spacetime',
  label: 'Spacetime diagram (2+1)',
  order: 4,
  description: 'Worldlines in a 2+1 spacetime diagram: the equatorial plane horizontal, coordinate time up. Light cones (45 degrees far away) narrow and tilt where the metric slows and drags light.',
  options: [
    { id: 'height', label: 'time axis', type: 'range', min: 0.3, max: 4, step: 0.1, value: 1.6, help: 'Height of the whole run relative to the view size' },
    { id: 'cones', label: 'light cones', type: 'checkbox', value: 1 },
    { id: 'coneSize', label: 'cone size', type: 'range', min: 0.01, max: 0.2, step: 0.005, value: 0.05 },
    { id: 'rings', label: 'cone rings', type: 'range', min: 1, max: 10, step: 1, value: 6 },
    { id: 'now', label: 'now plane', type: 'checkbox', value: 1 },
  ],
  available(data) {
    return data.trackSpecs.length ? null : 'needs at least one trajectory';
  },
  create(ctx) {
    const { THREE, data, options } = ctx;
    const extent = ctx.extent;
    const group = new THREE.Group();
    const H = options.height * extent;
    const duration = Math.max(ctx.tEnd - ctx.tStart, 1e-30);
    const yOf = (t) => (t - ctx.tStart) / duration * H;
    const notes = ['Cones: 45 degrees = light far from the mass, independent of the time-axis scale.'];

    // Time axis with ticks.
    {
      const pts = [new THREE.Vector3(0, 0, 0), new THREE.Vector3(0, H * 1.04, 0)];
      for (let k = 0; k <= 10; k++) {
        const y = H * k / 10;
        pts.push(new THREE.Vector3(-0.02 * extent, y, 0), new THREE.Vector3(0.02 * extent, y, 0));
      }
      const geom = new THREE.BufferGeometry().setFromPoints(pts);
      group.add(new THREE.LineSegments(geom, new THREE.LineBasicMaterial({ color: 0x8a93a6 })));
    }
    // Spheres become cylinders (their world tubes).
    for (const s of data.spheres) {
      const color = new THREE.Color(s.color || '#888888');
      const horizon = /horizon/i.test(s.name);
      const cyl = new THREE.Mesh(new THREE.CylinderGeometry(s.radius, s.radius, H, 64, 1, true),
        new THREE.MeshBasicMaterial({ color: horizon ? 0x000000 : color, transparent: true, opacity: horizon ? 0.75 : 0.18, side: THREE.DoubleSide, depthWrite: false }));
      cyl.position.y = H / 2;
      group.add(cyl);
      // Outline: rings at both ends and a few generators (an edges geometry would draw every facet).
      const pts = [];
      for (let m = 0; m < 96; m++) {
        const a0 = 2 * Math.PI * m / 96, a1 = 2 * Math.PI * (m + 1) / 96;
        for (const y of [0, H]) pts.push(new THREE.Vector3(s.radius * Math.cos(a0), y, s.radius * Math.sin(a0)), new THREE.Vector3(s.radius * Math.cos(a1), y, s.radius * Math.sin(a1)));
      }
      for (let m = 0; m < 8; m++) {
        const a = 2 * Math.PI * m / 8;
        pts.push(new THREE.Vector3(s.radius * Math.cos(a), 0, s.radius * Math.sin(a)), new THREE.Vector3(s.radius * Math.cos(a), H, s.radius * Math.sin(a)));
      }
      group.add(new THREE.LineSegments(new THREE.BufferGeometry().setFromPoints(pts),
        new THREE.LineBasicMaterial({ color: horizon ? 0x777777 : color, transparent: true, opacity: 0.45 })));
    }
    // Floor disk at t0.
    {
      const floor = new THREE.Mesh(new THREE.CircleGeometry(extent, 96), new THREE.MeshBasicMaterial({ color: 0x151a24, side: THREE.DoubleSide, transparent: true, opacity: 0.6 }));
      floor.rotation.x = -Math.PI / 2;
      group.add(floor);
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

    // Light cones on the now plane, from the slice's coordinate speeds of light.
    const coneGroup = new THREE.Group();
    const slice = data.slice;
    if (slice && slice.fields.coneOut) {
      const get = (name) => ctx.sampler(slice, name, 'coneValid');
      const out = get('coneOut'), inw = get('coneIn'), pro = get('conePro'), retro = get('coneRetro');
      const rDead = Math.max(slice.rMin || 0, 1e-6 * extent);
      // Innermost radius with valid cones along +x.
      let rValid = rDead;
      for (let R = rDead; R < extent; R += extent / 400) if (out(R, 0) !== null) { rValid = R; break; }
      const h = options.coneSize * extent;
      const rings = Math.round(options.rings);
      const rLow = Math.max(rValid * 1.08, 0.02 * extent), rHigh = 0.85 * extent;
      const pts = [], cols = [];
      let sMin = Infinity;
      const nodes = [];
      for (let k = 0; k < rings; k++) {
        const R = rings === 1 ? 0.5 * (rLow + rHigh) : rLow * Math.pow(rHigh / rLow, k / (rings - 1));
        const spokes = 8;
        for (let m = 0; m < spokes; m++) {
          const a = 2 * Math.PI * (m + 0.5 * (k % 2)) / spokes;
          const X = R * Math.cos(a), Y = R * Math.sin(a);
          const so = out(X, Y), si = inw(X, Y), sp = pro(X, Y), sr = retro(X, Y);
          if ([so, si, sp, sr].some(v => v === null)) continue;
          nodes.push({ X, Y, so, si, sp, sr });
          sMin = Math.min(sMin, so);
        }
      }
      const segs = 32;
      for (const nd of nodes) {
        const rx = nd.X / Math.hypot(nd.X, nd.Y), ry = nd.Y / Math.hypot(nd.X, nd.Y);
        const tx = -ry, ty = rx;
        const base = [];
        for (let m = 0; m <= segs; m++) {
          const psi = 2 * Math.PI * m / segs, c = Math.cos(psi), s = Math.sin(psi);
          const offR = c >= 0 ? nd.so * c : -nd.si * c;
          const offT = s >= 0 ? nd.sp * s : -nd.sr * s;
          base.push([nd.X + h * (offR * rx + offT * tx), nd.Y + h * (offR * ry + offT * ty)]);
        }
        const col = ctx.colormap(nd.so);
        const push = (a, ya, b, yb) => { pts.push(a[0], ya, -a[1], b[0], yb, -b[1]); cols.push(col.r, col.g, col.b, col.r, col.g, col.b); };
        for (let m = 0; m < segs; m++) push(base[m], h, base[m + 1], h);
        for (let m = 0; m < segs; m += 4) push([nd.X, nd.Y], 0, base[m], h);
      }
      const geom = new THREE.BufferGeometry();
      geom.setAttribute('position', new THREE.Float32BufferAttribute(pts, 3));
      geom.setAttribute('color', new THREE.Float32BufferAttribute(cols, 3));
      coneGroup.add(new THREE.LineSegments(geom, new THREE.LineBasicMaterial({ vertexColors: true })));
      ctx.drawScale(0, 1, 'outward light speed (cone colour)');
      if (Number.isFinite(sMin)) notes.push(`Slowest outward light shown: ${ctx.fmt(sMin)} c (coordinate speed).`);
    } else {
      notes.push('No slice in the scene: light cones need Visualization.writeSpacetimeSlice.');
    }
    coneGroup.visible = !!options.cones;
    group.add(coneGroup);
    ctx.setNotes(notes.join('\n'));

    return {
      group,
      place(x) {
        const p = ctx.planeXYZ(x);
        const v = ctx.toWorld(p.X, p.Y, 0);
        v.y = yOf(x[0] / (data.header.c || 1));
        return v;
      },
      update(clock) {
        const y = yOf(clock);
        nowGroup.position.y = y;
        coneGroup.position.y = y;
      },
      camera() {
        return { position: new THREE.Vector3(0.35 * extent, 0.85 * H, 2.1 * extent), target: new THREE.Vector3(0, 0.45 * H, 0) };
      },
      readout() {
        const lines = [`time axis: ${ctx.fmt(duration)} of coordinate time over ${ctx.fmt(H)} lengths`];
        return lines;
      },
      onOption(id, v) {
        if (id === 'cones') { coneGroup.visible = !!v; return true; }
        if (id === 'now') { nowGroup.visible = !!v; return true; }
        return false;
      },
    };
  },
});
