// Clocks style: a time-contracting grid. Every node of the equatorial slice carries a clock of an observer
// at rest in the slicing; it runs at alpha = d tau / dt (the lapse). The sheet sinks by the time each clock
// has lost against clocks far away, (1 - alpha)(t - t0), so it deepens as the simulation runs; clock faces
// tick at their local rate and each body carries a clock showing its own proper time.
(window.GRStyles = window.GRStyles || []).push({
  id: 'clocks',
  label: 'Clocks / time grid',
  order: 3,
  description: 'Gravitational time dilation: each grid point holds a clock running at the lapse alpha; the sheet sinks by the time it has lost against far-away clocks, (1 - alpha)(t - t0). Bodies carry their own proper-time clocks.',
  options: [
    { id: 'warp', label: 'sink depth', type: 'range', min: 0, max: 3, step: 0.05, value: 1, help: 'Depth of the deepest point at the end of the run, relative to the view size' },
    { id: 'clocks', label: 'clock faces', type: 'checkbox', value: 1 },
    { id: 'stride', label: 'clock spacing', type: 'range', min: 3, max: 20, step: 1, value: 8, help: 'Grid nodes between clock faces' },
    { id: 'turns', label: 'hand turns', type: 'range', min: 1, max: 40, step: 1, value: 8, help: 'Turns of a far-away clock hand over the whole run' },
    { id: 'grid', label: 'grid lines', type: 'checkbox', value: 1 },
  ],
  available(data) {
    return data.slice && data.slice.fields.timeLag ? null : 'needs an equatorial slice with the lapse (Visualization.writeSpacetimeSlice)';
  },
  create(ctx) {
    const { THREE, data, options } = ctx;
    const slice = data.slice;
    const extent = ctx.extent;
    const group = new THREE.Group();
    const nu = slice.u.length, nv = slice.v.length;
    const lagF = slice.fields.timeLag, validF = slice.fields.valid;
    const lagAt = ctx.sampler(slice, 'timeLag', 'valid');

    let lagMax = 0;
    for (let i = 0; i < nu; i++) for (let j = 0; j < nv; j++) {
      if (validF[i][j] > 0.5 && Math.hypot(slice.u[i], slice.v[j]) <= extent) lagMax = Math.max(lagMax, lagF[i][j]);
    }
    lagMax = lagMax || 1;
    // Weak fields (Earth: 1 - alpha ~ 7e-10) get exaggerated hands; the sheet is always normalized.
    const handBoost = lagMax < 0.05 ? 0.5 / lagMax : 1;
    const duration = Math.max(ctx.tEnd - ctx.tStart, 1e-30);
    const depthScale = () => options.warp * 0.45 * extent / lagMax;
    const progress = (t) => (t - ctx.tStart);                      // coordinate time elapsed
    const depthAt = (lag, t) => -depthScale() * lag * progress(t) / duration;
    // Colour by sqrt of the lag, so the gentle far field is not one flat colour.
    ctx.drawScale(0, lagMax, 'time lag 1 - alpha (sqrt colour)', (x) => ctx.colormap(1 - Math.sqrt(x)));
    const notes = ['Depth = time lost by a static clock since the start, (1 - alpha)(t - t0), scaled so the deepest point reaches the chosen depth at the end.'];
    if (handBoost > 1) notes.push(`Weak field: clock-hand slowdown exaggerated x${handBoost.toPrecision(3)} (largest 1 - alpha = ${lagMax.toExponential(2)}).`);
    ctx.setNotes(notes.join('\n'));

    // --- sheet -----------------------------------------------------------------------------------
    const positions = new Float32Array(nu * nv * 3), colors = new Float32Array(nu * nv * 3);
    const indices = [], gridIdx = [];
    const ok = (i, j) => validF[i][j] > 0.5 && Math.hypot(slice.u[i], slice.v[j]) <= extent * 1.0001;
    for (let i = 0; i < nu; i++) for (let j = 0; j < nv; j++) {
      const k = i * nv + j;
      positions[3 * k] = slice.u[i]; positions[3 * k + 2] = -slice.v[j];
      const c = ok(i, j) ? ctx.colormap(1 - Math.sqrt(lagF[i][j] / lagMax)) : new THREE.Color(0, 0, 0);
      colors[3 * k] = c.r * 0.8; colors[3 * k + 1] = c.g * 0.8; colors[3 * k + 2] = c.b * 0.8;
    }
    for (let i = 0; i < nu - 1; i++) for (let j = 0; j < nv - 1; j++) {
      if (!(ok(i, j) && ok(i + 1, j) && ok(i, j + 1) && ok(i + 1, j + 1))) continue;
      const a = i * nv + j, b = a + 1, c = a + nv, d = c + 1;
      indices.push(a, b, c, b, d, c);
    }
    for (let i = 0; i < nu; i += 2) for (let j = 0; j < nv - 1; j++) if (ok(i, j) && ok(i, j + 1)) gridIdx.push(i * nv + j, i * nv + j + 1);
    for (let j = 0; j < nv; j += 2) for (let i = 0; i < nu - 1; i++) if (ok(i, j) && ok(i + 1, j)) gridIdx.push(i * nv + j, (i + 1) * nv + j);
    const posAttr = new THREE.BufferAttribute(positions, 3);
    const sheetGeom = new THREE.BufferGeometry();
    sheetGeom.setAttribute('position', posAttr);
    sheetGeom.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    sheetGeom.setIndex(indices);
    const sheet = new THREE.Mesh(sheetGeom, new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide, transparent: true, opacity: 0.92 }));
    group.add(sheet);
    const gridGeom = new THREE.BufferGeometry();
    gridGeom.setAttribute('position', posAttr);
    gridGeom.setIndex(gridIdx);
    const grid = new THREE.LineSegments(gridGeom, new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.22 }));
    grid.visible = !!options.grid;
    group.add(grid);

    // --- clock faces ---------------------------------------------------------------------------------
    const stride = Math.round(options.stride);
    const spacing = (slice.u[nu - 1] - slice.u[0]) / (nu - 1);
    const faces = [];
    for (let i = Math.floor(stride / 2); i < nu; i += stride) for (let j = Math.floor(stride / 2); j < nv; j += stride) {
      if (ok(i, j)) faces.push({ i, j, X: slice.u[i], Y: slice.v[j], lag: lagF[i][j] });
    }
    const faceR = 0.36 * stride * spacing, ringSegs = 28;
    const ringPos = new Float32Array(faces.length * ringSegs * 6), handPos = new Float32Array(faces.length * 6 + 64 * 6);
    const ringCol = new Float32Array(faces.length * ringSegs * 6);
    faces.forEach((f, q) => {
      const c = ctx.colormap(1 - Math.sqrt(f.lag / lagMax));
      for (let m = 0; m < ringSegs; m++) for (let e = 0; e < 2; e++) {
        const o = (q * ringSegs + m) * 6 + 3 * e;
        ringCol[o] = 0.4 + 0.6 * c.r; ringCol[o + 1] = 0.4 + 0.6 * c.g; ringCol[o + 2] = 0.4 + 0.6 * c.b;
      }
    });
    const ringGeom = new THREE.BufferGeometry();
    ringGeom.setAttribute('position', new THREE.BufferAttribute(ringPos, 3));
    ringGeom.setAttribute('color', new THREE.BufferAttribute(ringCol, 3));
    const rings = new THREE.LineSegments(ringGeom, new THREE.LineBasicMaterial({ vertexColors: true }));
    const handGeom = new THREE.BufferGeometry();
    handGeom.setAttribute('position', new THREE.BufferAttribute(handPos, 3));
    const hands = new THREE.LineSegments(handGeom, new THREE.LineBasicMaterial({ color: 0xffffff }));
    const faceGroup = new THREE.Group();
    faceGroup.add(rings, hands);
    faceGroup.visible = !!options.clocks;
    group.add(faceGroup);

    // Bodies' own clocks: one ring and hand each, drawn around the body.
    const bodyFaces = new THREE.Group();
    group.add(bodyFaces);
    const bodyRingSegs = 32;
    const bodyClocks = ctx.data.trackSpecs.map((spec) => {
      const color = new THREE.Color(spec.color || '#ffffff');
      const ring = new THREE.LineLoop(new THREE.BufferGeometry().setFromPoints(
        Array.from({ length: bodyRingSegs }, (_, m) => new THREE.Vector3(Math.cos(2 * Math.PI * m / bodyRingSegs), 0, Math.sin(2 * Math.PI * m / bodyRingSegs)))),
        new THREE.LineBasicMaterial({ color }));
      const hand = new THREE.Line(new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(), new THREE.Vector3(1, 0, 0)]), new THREE.LineBasicMaterial({ color: 0xffffff }));
      const g = new THREE.Group();
      g.add(ring, hand);
      bodyFaces.add(g);
      return { group: g, hand, lambda0: spec.samples[0].lambda };
    });
    const bodyFaceR = 1.4 * faceR;

    const handAngle = (elapsedOwn, t) => {
      // Exaggerate only the lag against far-away clocks, so far clocks keep their true rate.
      const far = progress(t);
      const shown = far - handBoost * (far - elapsedOwn);
      return 2 * Math.PI * options.turns * shown / duration;
    };

    function refresh(t) {
      for (let i = 0; i < nu; i++) for (let j = 0; j < nv; j++) {
        const k = i * nv + j;
        positions[3 * k + 1] = ok(i, j) ? depthAt(lagF[i][j], t) : depthAt(lagMax, t);
      }
      posAttr.needsUpdate = true;
      if (options.clocks) {
        const lift = 0.004 * extent;
        faces.forEach((f, q) => {
          const y = depthAt(f.lag, t) + lift;
          for (let m = 0; m < ringSegs; m++) {
            const a0 = 2 * Math.PI * m / ringSegs, a1 = 2 * Math.PI * (m + 1) / ringSegs, o = (q * ringSegs + m) * 6;
            ringPos[o] = f.X + faceR * Math.cos(a0); ringPos[o + 1] = y; ringPos[o + 2] = -(f.Y + faceR * Math.sin(a0));
            ringPos[o + 3] = f.X + faceR * Math.cos(a1); ringPos[o + 4] = y; ringPos[o + 5] = -(f.Y + faceR * Math.sin(a1));
          }
          // A static clock's own time is alpha (t - t0).
          const ang = handAngle((1 - f.lag) * progress(t), t);
          const o = q * 6;
          handPos[o] = f.X; handPos[o + 1] = y; handPos[o + 2] = -f.Y;
          handPos[o + 3] = f.X + 0.85 * faceR * Math.sin(ang); handPos[o + 4] = y; handPos[o + 5] = -(f.Y + 0.85 * faceR * Math.cos(ang));
        });
        ringGeom.attributes.position.needsUpdate = true;
        handGeom.attributes.position.needsUpdate = true;
        handGeom.setDrawRange(0, faces.length * 2);
      }
      const bodies = ctx.bodiesNow();
      bodyClocks.forEach((bc, q) => {
        const b = bodies[q];
        if (!b || !b.position) return;
        bc.group.position.copy(b.position);
        bc.group.scale.setScalar(bodyFaceR);
        const ang = handAngle(b.sample.lambda - bc.lambda0, t);
        const attr = bc.hand.geometry.getAttribute('position');
        attr.setXYZ(1, Math.sin(ang), 0, -Math.cos(ang));
        attr.needsUpdate = true;
      });
    }

    const sheetHeight = (X, Y, t) => {
      const lag = lagAt(X, Y);
      return depthAt(lag === null ? lagMax : lag, t);
    };

    return {
      group,
      dynamicPlacement: true,
      place(x, sample, t) {
        const p = ctx.planeXYZ(x);
        return ctx.toWorld(p.X, p.Y, sheetHeight(p.X, p.Y, t) + 0.012 * extent);
      },
      update(clock) { refresh(clock); },
      camera() {
        return { position: ctx.toWorld(0, -1.35 * extent, 0.75 * extent), target: new THREE.Vector3(0, -0.15 * extent * options.warp, 0) };
      },
      readout() {
        const t = ctx.clock, lines = [];
        lines.push(`deepest clock: 1 - alpha = ${ctx.fmt(lagMax)}, has lost ${ctx.fmt(lagMax * progress(t))}`);
        for (const b of ctx.bodiesNow()) {
          const own = b.sample.lambda - (bodyClocks[ctx.bodiesNow().indexOf(b)] || { lambda0: 0 }).lambda0;
          lines.push(`${b.name}: own time ${ctx.fmt(own)}, behind far clocks by ${ctx.fmt(progress(t) - own)}`);
        }
        return lines;
      },
      onOption(id, v) {
        if (id === 'warp' || id === 'turns') { refresh(ctx.clock); return true; }
        if (id === 'grid') { grid.visible = !!v; return true; }
        if (id === 'clocks') { faceGroup.visible = !!v; refresh(ctx.clock); return true; }
        return false;
      },
    };
  },
});
