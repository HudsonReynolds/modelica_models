// River style (Hamilton & Lisle, Am. J. Phys. 76, 519 (2008)): space flows into the hole at the speed of
// free fall from rest at infinity, v = sqrt(1 + g_tt) (sqrt(rs/r) for Schwarzschild), swirling with the
// frame-dragging rate for a spinning hole. Light moves at c relative to the flowing space, so a flash
// emitted inside the horizon is carried inward. Drawn in the equatorial plane: one spatial dimension is
// dropped. Gullstrand-Painleve slices are flat, so this flat picture is exact for Schwarzschild.
(window.GRStyles = window.GRStyles || []).push({
  id: 'river',
  label: 'River model',
  order: 2,
  description: 'Space flows inward at the speed of free fall from infinity (river model). Streaks move with the flowing space; light flashes expand at c relative to it, so inside the horizon even light is swept in. Equatorial plane only.',
  options: [
    { id: 'tracers', label: 'tracers', type: 'range', min: 200, max: 6000, step: 100, value: 2500 },
    { id: 'flow', label: 'flow rate', type: 'range', min: 0.1, max: 4, step: 0.1, value: 1, help: 'Animation rate of the river; at 1, light crosses 1/8 of the view per second' },
    { id: 'flashes', label: 'light flashes', type: 'checkbox', value: 1 },
    { id: 'rings', label: 'flash rings', type: 'range', min: 1, max: 10, step: 1, value: 5, help: 'Rings of flash sources' },
    { id: 'map', label: 'speed map', type: 'checkbox', value: 1 },
  ],
  available(data) {
    return data.slice && data.slice.fields.river ? null : 'needs an equatorial slice (Visualization.writeSpacetimeSlice, or SceneRecorder with writeSlice)';
  },
  create(ctx) {
    const { THREE, data, options } = ctx;
    const slice = data.slice;
    const extent = ctx.extent;
    const group = new THREE.Group();
    const river = ctx.sampler(slice, 'river');
    const drag = slice.fields.drag ? ctx.sampler(slice, 'drag') : () => 0;
    const rDead = Math.max(slice.rMin || 0, 1e-6 * extent);

    // Speeds over the evaluated part of the grid.
    let vMax = 0, vOuter = 0;
    slice.u.forEach((u, i) => slice.v.forEach((v, j) => {
      const R = Math.hypot(u, v);
      if (R <= rDead || R > extent) return;
      vMax = Math.max(vMax, slice.fields.river[i][j]);
    }));
    vOuter = vMax;
    // Weak fields (Earth: 4e-5 c) are exaggerated so the flow is visible; strong fields are drawn true to c.
    const boost = vOuter > 0 && vOuter < 0.2 ? 0.25 / vOuter : 1;
    const colourTop = boost > 1 ? vOuter : Math.min(Math.max(vMax, 1.2), 2.0);
    ctx.drawScale(0, colourTop, 'river speed / c');
    const notes = [];
    if (boost > 1) notes.push(`Weak field: river speeds exaggerated x${boost.toPrecision(3)} (fastest ${vOuter.toExponential(2)} c); light flashes are not.`);
    if (vMax > 1) notes.push('White contour: the river reaches c (horizon; static limit for a spinning hole).');
    ctx.setNotes(notes.join('\n'));

    // Velocity of the flowing space at (X, Y) in display length per unit ct, or null where undefined.
    const velocity = (X, Y, out) => {
      const R = Math.hypot(X, Y);
      if (R <= rDead * 1.0001 || R >= extent) return null;
      const v = river(X, Y);
      if (v === null) return null;
      const d = drag(X, Y) || 0;
      out.x = boost * (-v * X / R - d * Y / R);
      out.y = boost * (-v * Y / R + d * X / R);
      return out;
    };

    // --- background: speed map ---------------------------------------------------------------------
    {
      const nu = slice.u.length, nv = slice.v.length;
      const positions = [], colors = [], indices = [];
      for (let i = 0; i < nu; i++) for (let j = 0; j < nv; j++) {
        const X = slice.u[i], Y = slice.v[j], R = Math.hypot(X, Y);
        positions.push(X, -0.002 * extent, -Y);
        const dead = R <= rDead;
        const c = dead ? new THREE.Color(0, 0, 0) : ctx.colormap(slice.fields.river[i][j] / colourTop).multiplyScalar(0.55);
        const fade = R > extent ? 0.25 : 1;
        colors.push(c.r * fade, c.g * fade, c.b * fade);
      }
      for (let i = 0; i < nu - 1; i++) for (let j = 0; j < nv - 1; j++) {
        // Only cells inside the circle of radius extent, where tracers flow.
        if (Math.hypot(0.5 * (slice.u[i] + slice.u[i + 1]), 0.5 * (slice.v[j] + slice.v[j + 1])) > extent) continue;
        // Skip cells touching unevaluated nodes (inside a body or the innermost radius).
        if ([[i, j], [i + 1, j], [i, j + 1], [i + 1, j + 1]].some(([a, b]) => Math.hypot(slice.u[a], slice.v[b]) <= rDead)) continue;
        const a = i * nv + j, b = a + 1, c = a + nv, d = c + 1;
        indices.push(a, b, c, b, d, c);
      }
      const geom = new THREE.BufferGeometry();
      geom.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
      geom.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
      geom.setIndex(indices);
      const map = new THREE.Mesh(geom, new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide }));
      map.visible = !!options.map;
      map.name = 'map';
      group.add(map);
    }

    // --- contour where the river reaches c (marching squares) ---------------------------------------
    if (vMax > 1) {
      const f = slice.fields.river, u = slice.u, v = slice.v, pts = [];
      const lerp = (a, b, fa, fb) => a + (b - a) * (1 - fa) / (fb - fa);
      for (let i = 0; i < u.length - 1; i++) for (let j = 0; j < v.length - 1; j++) {
        if (Math.hypot(u[i], v[j]) <= rDead * 1.5) continue;
        const c00 = f[i][j], c10 = f[i + 1][j], c01 = f[i][j + 1], c11 = f[i + 1][j + 1];
        const cross = [];
        if ((c00 - 1) * (c10 - 1) < 0) cross.push([lerp(u[i], u[i + 1], c00, c10), v[j]]);
        if ((c10 - 1) * (c11 - 1) < 0) cross.push([u[i + 1], lerp(v[j], v[j + 1], c10, c11)]);
        if ((c01 - 1) * (c11 - 1) < 0) cross.push([lerp(u[i], u[i + 1], c01, c11), v[j + 1]]);
        if ((c00 - 1) * (c01 - 1) < 0) cross.push([u[i], lerp(v[j], v[j + 1], c00, c01)]);
        for (let k = 0; k + 1 < cross.length; k += 2) {
          pts.push(cross[k][0], 0.003 * extent, -cross[k][1], cross[k + 1][0], 0.003 * extent, -cross[k + 1][1]);
        }
      }
      const geom = new THREE.BufferGeometry();
      geom.setAttribute('position', new THREE.Float32BufferAttribute(pts, 3));
      group.add(new THREE.LineSegments(geom, new THREE.LineBasicMaterial({ color: 0xffffff })));
    }

    // --- reference spheres as circles; the horizon as a black disk ---------------------------------
    for (const s of data.spheres) {
      const ring = [];
      for (let k = 0; k <= 96; k++) {
        const a = 2 * Math.PI * k / 96;
        ring.push(new THREE.Vector3(s.radius * Math.cos(a), 0.004 * extent, -s.radius * Math.sin(a)));
      }
      const color = new THREE.Color(s.color || '#888888');
      group.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(ring),
        new THREE.LineBasicMaterial({ color: color.getHex() === 0 ? 0x666666 : color })));
      {
        const horizon = /horizon/i.test(s.name);
        const disk = new THREE.Mesh(new THREE.CircleGeometry(s.radius, 64),
          new THREE.MeshBasicMaterial({ color: horizon ? 0x000000 : color, transparent: true, opacity: horizon ? 0.65 : 0.35 }));
        disk.rotation.x = -Math.PI / 2;
        disk.position.y = 0.001 * extent;
        group.add(disk);
      }
    }

    // --- tracers --------------------------------------------------------------------------------
    const n = Math.round(options.tracers);
    const px = new Float32Array(n), py = new Float32Array(n);
    const spawn = (k, anywhere) => {
      const a = Math.random() * 2 * Math.PI;
      const rIn = rDead * 1.05;
      const R = anywhere ? Math.sqrt(rIn * rIn + Math.random() * (extent * extent - rIn * rIn)) : extent * (0.82 + 0.17 * Math.random());
      px[k] = R * Math.cos(a);
      py[k] = R * Math.sin(a);
    };
    for (let k = 0; k < n; k++) spawn(k, true);
    const streakPos = new Float32Array(n * 6), streakCol = new Float32Array(n * 6);
    const streakGeom = new THREE.BufferGeometry();
    streakGeom.setAttribute('position', new THREE.BufferAttribute(streakPos, 3));
    streakGeom.setAttribute('color', new THREE.BufferAttribute(streakCol, 3));
    const streaks = new THREE.LineSegments(streakGeom, new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, blending: THREE.AdditiveBlending, depthWrite: false }));
    group.add(streaks);

    // --- light flashes: circles of radius c * age whose centres drift with the river ---------------
    const ringsCount = Math.round(options.rings);
    const sources = [];
    const rLow = Math.max(rDead * 1.3, 0.03 * extent), rHigh = 0.85 * extent;
    for (let k = 0; k < ringsCount; k++) {
      const R = ringsCount === 1 ? 0.5 * (rLow + rHigh) : rLow * Math.pow(rHigh / rLow, k / (ringsCount - 1));
      const spokes = 8;
      for (let m = 0; m < spokes; m++) {
        const a = 2 * Math.PI * (m + 0.5 * (k % 2)) / spokes;
        sources.push({ X: R * Math.cos(a), Y: R * Math.sin(a), age: 0, cx: 0, cy: 0, alive: false, phase: Math.random() });
      }
    }
    const period = 0.12 * extent;            // flash lifetime in units of ct
    const segs = 40;
    const flashPos = new Float32Array(sources.length * segs * 6), flashCol = new Float32Array(sources.length * segs * 6);
    const flashGeom = new THREE.BufferGeometry();
    flashGeom.setAttribute('position', new THREE.BufferAttribute(flashPos, 3));
    flashGeom.setAttribute('color', new THREE.BufferAttribute(flashCol, 3));
    const flashes = new THREE.LineSegments(flashGeom, new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, blending: THREE.AdditiveBlending, depthWrite: false }));
    flashes.visible = !!options.flashes;
    group.add(flashes);
    let flowTime = 0;
    sources.forEach(s => { s.age = s.phase * period; s.cx = s.X; s.cy = s.Y; s.alive = true; });

    const tmp = { x: 0, y: 0 }, tmp2 = { x: 0, y: 0 };
    const advect = (x, y, h) => {
      // Midpoint rule; returns null when the point leaves the defined region.
      if (!velocity(x, y, tmp)) return null;
      const mx = x + 0.5 * h * tmp.x, my = y + 0.5 * h * tmp.y;
      if (!velocity(mx, my, tmp2)) return null;
      return [x + h * tmp2.x, y + h * tmp2.y];
    };

    function step(dt) {
      const h = dt * options.flow * extent / 8;   // ct advanced this frame
      flowTime += h;
      const lift = 0.006 * extent;
      for (let k = 0; k < n; k++) {
        let moved = h > 0 ? advect(px[k], py[k], h) : [px[k], py[k]];
        if (!moved) { spawn(k, false); moved = [px[k], py[k]]; }
        px[k] = moved[0]; py[k] = moved[1];
        const vel = velocity(px[k], py[k], tmp);
        const speed = vel ? Math.hypot(vel.x, vel.y) / boost : 0;
        // Streak = the displacement over the last 0.3 s of animation, capped.
        const ux = vel && speed > 0 ? vel.x / (speed * boost) : 0, uy = vel && speed > 0 ? vel.y / (speed * boost) : 0;
        const tail = Math.min(0.08 * extent, speed * boost * 0.3 * options.flow * extent / 8);
        const o = k * 6;
        streakPos[o] = px[k] - ux * tail; streakPos[o + 1] = lift; streakPos[o + 2] = -(py[k] - uy * tail);
        streakPos[o + 3] = px[k]; streakPos[o + 4] = lift; streakPos[o + 5] = -py[k];
        const c = ctx.colormap(0.35 + 0.65 * speed / colourTop);
        const hot = speed > 1 ? 1 : 0;
        streakCol[o] = 0.05; streakCol[o + 1] = 0.05; streakCol[o + 2] = 0.08;
        streakCol[o + 3] = hot ? 1.0 : c.r; streakCol[o + 4] = hot ? 0.45 : c.g; streakCol[o + 5] = hot ? 0.35 : c.b;
      }
      streakGeom.attributes.position.needsUpdate = true;
      streakGeom.attributes.color.needsUpdate = true;

      if (options.flashes) {
        let o = 0;
        for (const s of sources) {
          s.age += h;
          if (s.age >= period) { s.age = s.age % period; s.cx = s.X; s.cy = s.Y; s.alive = true; }
          else if (s.alive && h > 0) {
            const moved = advect(s.cx, s.cy, h);
            if (moved) { s.cx = moved[0]; s.cy = moved[1]; } else s.alive = false;
          }
          const radius = s.age;   // light: c = 1 in display length per ct
          const fade = s.alive ? Math.pow(1 - s.age / period, 1.5) : 0;
          for (let m = 0; m < segs; m++) {
            const a0 = 2 * Math.PI * m / segs, a1 = 2 * Math.PI * (m + 1) / segs;
            flashPos[o] = s.cx + radius * Math.cos(a0); flashPos[o + 1] = lift * 1.2; flashPos[o + 2] = -(s.cy + radius * Math.sin(a0));
            flashPos[o + 3] = s.cx + radius * Math.cos(a1); flashPos[o + 4] = lift * 1.2; flashPos[o + 5] = -(s.cy + radius * Math.sin(a1));
            for (let q = 0; q < 6; q += 3) { flashCol[o + q] = 0.55 * fade; flashCol[o + q + 1] = 0.85 * fade; flashCol[o + q + 2] = 1.0 * fade; }
            o += 6;
          }
        }
        flashGeom.attributes.position.needsUpdate = true;
        flashGeom.attributes.color.needsUpdate = true;
      }
    }
    step(0);

    return {
      group,
      place(x) {
        const p = ctx.planeXYZ(x);
        return ctx.toWorld(p.X, p.Y, 0.012 * extent);
      },
      update(clock, dt) { step(dt); },
      camera() {
        return { position: ctx.toWorld(0, -0.75 * extent, 1.7 * extent), target: new THREE.Vector3(0, 0, 0) };
      },
      readout() {
        const lines = [];
        for (const b of ctx.bodiesNow()) {
          const p = ctx.planeXYZ(b.sample.x);
          const v = river(p.X, p.Y);
          if (v !== null) lines.push(`${b.name}: space flows in at ${ctx.fmt(v)} c`);
        }
        return lines;
      },
      onOption(id, v) {
        if (id === 'flow') return true;
        if (id === 'flashes') { flashes.visible = !!v; return true; }
        if (id === 'map') { group.children.forEach(o => { if (o.name === 'map') o.visible = !!v; }); return true; }
        return false;
      },
    };
  },
});
