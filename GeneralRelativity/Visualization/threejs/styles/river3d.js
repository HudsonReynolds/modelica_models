// River 3D style: the river of space as a connected 3D grid falling into the hole. Free-fall parcels that
// start at the same radius stay on a common sphere (exact in spherical symmetry), so the grid is a stack of
// concentric shells, each a latitude/longitude net, joined by radial spokes. Every shell follows the same
// trajectory r(tau) with dr/dtau = -v(r), v = sqrt(1 + g_tt), offset in phase, so shells stream inward
// continuously and are re-emitted at the outer edge; frame dragging (Kerr) twists each shell about the
// spin axis at omega(r) = drag / r. Speeds come from the equatorial slice, averaged over azimuth.
(window.GRStyles = window.GRStyles || []).push({
  id: 'river3d',
  label: 'River model 3D',
  order: 2.5,
  description: 'The river of space as a 3D grid: shells of free-falling space stream inward, stretching radially as the flow speeds up, crossing the horizon at the speed of light and twisting with frame dragging. A wedge is cut away to show the inside.',
  options: [
    { id: 'shells', label: 'shells', type: 'range', min: 4, max: 40, step: 1, value: 10 },
    { id: 'meridians', label: 'meridians', type: 'range', min: 4, max: 36, step: 2, value: 12 },
    { id: 'parallels', label: 'parallels', type: 'range', min: 3, max: 18, step: 1, value: 5 },
    { id: 'spokes', label: 'radial lines', type: 'checkbox', value: 1 },
    { id: 'cut', label: 'cut-away', type: 'select', choices: [{ value: 0, label: 'none' }, { value: 1, label: 'front wedge' }, { value: 2, label: 'upper half' }], value: 1 },
    { id: 'outer', label: 'outer radius', type: 'range', min: 0.2, max: 1, step: 0.05, value: 0.7, help: 'Fraction of the scene size' },
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
    // Stay 1.5 grid cells outside the unevaluated core: bilinear samples there would mix in its zeros.
    const rDead = Math.max(slice.rMin || 0, 1e-6 * extent) + 1.5 * (slice.u[slice.u.length - 1] - slice.u[0]) / (slice.u.length - 1);
    // At least a little outside the core, whatever the slider says (a planet can fill much of the view).
    const rOuter = Math.max(options.outer * extent, 1.3 * rDead);

    // Azimuthally averaged equatorial profiles v(r) and omega(r).
    const nProfile = 400, rP = [], vP = [], wP = [];
    for (let k = 0; k < nProfile; k++) {
      const r = rDead * 1.02 + (rOuter - rDead * 1.02) * k / (nProfile - 1);
      let v = 0, w = 0, n = 0;
      for (let m = 0; m < 16; m++) {
        const a = 2 * Math.PI * m / 16, X = r * Math.cos(a), Y = r * Math.sin(a);
        const vv = river(X, Y), dd = drag(X, Y);
        if (vv === null) continue;
        v += vv; w += (dd || 0) / r; n++;
      }
      if (n) { rP.push(r); vP.push(v / n); wP.push(w / n); }
    }
    const lookup = (arr, r) => {
      if (r <= rP[0]) return arr[0];
      if (r >= rP[rP.length - 1]) return arr[arr.length - 1];
      const f = (r - rP[0]) / (rP[rP.length - 1] - rP[0]) * (rP.length - 1), i = Math.floor(f);
      return arr[i] + (f - i) * (arr[i + 1] - arr[i]);
    };
    const vMax = Math.max(...vP), vMin = Math.min(...vP);
    const colourTop = vMax > 0.2 ? Math.min(Math.max(vMax, 1.2), 2) : vMax;

    // One fall trajectory from rOuter to rDead, in flow time s (units of ct), with the accumulated twist.
    const traj = { s: [0], r: [rOuter], psi: [0] };
    {
      let r = rOuter, s = 0, psi = 0;
      const vFloor = Math.max(vMin, 1e-12);
      for (let k = 0; k < 20000 && r > rDead * 1.03; k++) {
        const v = Math.max(lookup(vP, r), vFloor);
        const ds = 0.004 * Math.max(r - rDead, 0.02 * extent) / v;
        const rm = r - 0.5 * ds * v;
        const vm = Math.max(lookup(vP, rm), vFloor);
        r -= ds * vm;
        psi += ds * lookup(wP, rm);
        s += ds;
        traj.s.push(s); traj.r.push(r); traj.psi.push(psi);
      }
    }
    const period = traj.s[traj.s.length - 1] || 1;
    const along = (s) => {
      s = ((s % period) + period) % period;
      let lo = 0, hi = traj.s.length - 1;
      while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (traj.s[mid] <= s) lo = mid; else hi = mid; }
      const f = (s - traj.s[lo]) / Math.max(traj.s[hi] - traj.s[lo], 1e-30);
      return { r: traj.r[lo] + f * (traj.r[hi] - traj.r[lo]), psi: traj.psi[lo] + f * (traj.psi[hi] - traj.psi[lo]) };
    };

    const nShell = Math.round(options.shells), nMer = Math.round(options.meridians), nPar = Math.round(options.parallels);
    const segPar = 48, segMer = 24;
    // Unit directions: parallels (theta fixed) and meridians (phi fixed); twist rotates about z.
    const thetaPar = Array.from({ length: nPar }, (_, k) => Math.PI * (k + 1) / (nPar + 1));
    const phiMer = Array.from({ length: nMer }, (_, k) => 2 * Math.PI * k / nMer);
    const cutAway = (theta, phi) => {
      if (options.cut === 1) {
        // 120-degree wedge facing the default camera (which looks from -Y).
        const p = Math.atan2(Math.sin(phi), Math.cos(phi));
        return p > -5 * Math.PI / 6 && p < -Math.PI / 6;
      }
      if (options.cut === 2) return theta < Math.PI / 2 - 1e-9;
      return false;
    };
    const maxSegs = nShell * (nPar * segPar + nMer * segMer) + nShell * nPar * nMer;
    const pos = new Float32Array(maxSegs * 6), col = new Float32Array(maxSegs * 6);
    const geom = new THREE.BufferGeometry();
    geom.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    geom.setAttribute('color', new THREE.BufferAttribute(col, 3));
    const lines = new THREE.LineSegments(geom, new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.9 }));
    group.add(lines);

    // Spin axis and the speed-of-light sphere (where the averaged river reaches c).
    group.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(0, -rOuter, 0), new THREE.Vector3(0, rOuter, 0)]),
      new THREE.LineBasicMaterial({ color: 0x555b68 })));
    let rLight = null;
    // Outermost radius where the river reaches c (searched from outside in).
    for (let k = rP.length - 1; k > 0; k--) if ((vP[k - 1] - 1) * (vP[k] - 1) <= 0 && vP[k - 1] !== vP[k]) { rLight = rP[k - 1] + (1 - vP[k - 1]) / (vP[k] - vP[k - 1]) * (rP[k] - rP[k - 1]); break; }
    for (const s of data.spheres) {
      const horizon = /horizon/i.test(s.name);
      group.add(new THREE.Mesh(new THREE.SphereGeometry(s.radius, 48, 32),
        new THREE.MeshBasicMaterial({ color: horizon ? 0x000000 : new THREE.Color(s.color || '#888888'), transparent: !horizon, opacity: horizon ? 1 : 0.5 })));
    }
    if (rLight) {
      group.add(new THREE.Mesh(new THREE.SphereGeometry(rLight, 48, 32),
        new THREE.MeshBasicMaterial({ color: 0xff6b5a, wireframe: true, transparent: true, opacity: 0.18 })));
    }
    ctx.drawScale(0, colourTop, 'river speed / c (shell colour)');
    const clockT = ctx.freeFallClock();
    const cLight = data.header.c || 1;
    const notes = [`A shell takes ${ctx.fmt(period)} of free-fall time (units of ct) to fall from r = ${ctx.fmt(rOuter)} to ${ctx.fmt(rDead)}. The shells move on the simulation clock, in free-fall time, so a body dropped from rest far away falls with them.`];
    if (rLight) notes.push(`Red wire sphere: the river reaches c at r = ${ctx.fmt(rLight)}; shells keep falling, faster than light, inside it.`);
    if (vMax < 0.2) notes.push(`Weak field: fastest flow ${vMax.toExponential(2)} c. Only the timing is scaled; the shell spacing is true to v(r).`);
    if (wP.some(w => Math.abs(w) > 0)) notes.push('Shells twist about the spin axis with the frame-dragging rate (equatorial value at each radius).');
    ctx.setNotes(notes.join('\n'));

    let flow = 0;
    const colourOf = (r) => {
      const v = lookup(vP, r);
      if (v > 1) return new THREE.Color(1.0, 0.42 + 0.3 * Math.max(0, 2 - v), 0.35);
      return ctx.colormap(0.25 + 0.75 * v / colourTop);
    };
    function rebuild() {
      const shells = [];
      for (let k = 0; k < nShell; k++) shells.push(along(flow + period * k / nShell));
      shells.sort((a, b) => b.r - a.r);
      let o = 0;
      const put = (x0, y0, z0, x1, y1, z1, c0, c1) => {
        // Picture (x, y, z) with z the spin axis -> world (x, z, -y).
        pos[o] = x0; pos[o + 1] = z0; pos[o + 2] = -y0; pos[o + 3] = x1; pos[o + 4] = z1; pos[o + 5] = -y1;
        col[o] = c0.r; col[o + 1] = c0.g; col[o + 2] = c0.b; col[o + 3] = c1.r; col[o + 4] = c1.g; col[o + 5] = c1.b;
        o += 6;
      };
      const point = (r, theta, phi) => [r * Math.sin(theta) * Math.cos(phi), r * Math.sin(theta) * Math.sin(phi), r * Math.cos(theta)];
      for (const sh of shells) {
        const c = colourOf(sh.r);
        // Fade shells in at the outer edge and out at the centre so the recycling is not visible.
        const fade = Math.min(1, (rOuter - sh.r) / (0.08 * rOuter) + 0.05, (sh.r - rDead) / (0.15 * Math.max(rDead, 0.02 * extent)));
        const cf = c.clone().multiplyScalar(Math.max(0, fade));
        for (const th of thetaPar) for (let m = 0; m < segPar; m++) {
          const p0 = 2 * Math.PI * m / segPar + sh.psi, p1 = 2 * Math.PI * (m + 1) / segPar + sh.psi;
          if (cutAway(th, 0.5 * (p0 + p1))) continue;
          const a = point(sh.r, th, p0), b = point(sh.r, th, p1);
          put(a[0], a[1], a[2], b[0], b[1], b[2], cf, cf);
        }
        for (const ph of phiMer) {
          if (cutAway(Math.PI / 2, ph + sh.psi) && options.cut === 1) continue;
          for (let m = 0; m < segMer; m++) {
            const t0 = Math.PI * m / segMer, t1 = Math.PI * (m + 1) / segMer;
            if (cutAway(0.5 * (t0 + t1), ph + sh.psi)) continue;
            const a = point(sh.r, t0, ph + sh.psi), b = point(sh.r, t1, ph + sh.psi);
            put(a[0], a[1], a[2], b[0], b[1], b[2], cf, cf);
          }
        }
      }
      if (options.spokes) {
        // Radial lines join the same (theta, phi) parcel on consecutive shells.
        for (let k = 0; k + 1 < shells.length; k++) {
          const A = shells[k], B = shells[k + 1];
          if (A.r < B.r) continue;
          const fa = Math.min(1, (rOuter - A.r) / (0.08 * rOuter) + 0.05), fb = Math.min(1, (B.r - rDead) / (0.15 * Math.max(rDead, 0.02 * extent)));
          const ca = colourOf(A.r).multiplyScalar(Math.max(0, fa) * 0.7), cb = colourOf(B.r).multiplyScalar(Math.max(0, Math.min(fa, fb)) * 0.7);
          for (const th of thetaPar) for (const ph of phiMer) {
            if (cutAway(th, ph + A.psi)) continue;
            const a = point(A.r, th, ph + A.psi), b = point(B.r, th, ph + B.psi);
            put(a[0], a[1], a[2], b[0], b[1], b[2], ca, cb);
          }
        }
      }
      geom.setDrawRange(0, o / 3);
      geom.attributes.position.needsUpdate = true;
      geom.attributes.color.needsUpdate = true;
    }
    rebuild();

    return {
      group,
      place(x) {
        const p = ctx.planeXYZ(x);
        return ctx.toWorld(p.X, p.Y, p.Z);
      },
      timeOf: clockT ? clockT.timeOf : undefined,
      clockLabel: clockT ? 'T' : 't',
      update(clock) {
        // ct of the simulation clock; the trajectory is periodic in it, so jumps (scrubbing, looping) are exact.
        flow = cLight * clock;
        rebuild();
      },
      camera() {
        return { position: ctx.toWorld(0.35 * rOuter, -1.9 * rOuter, 1.0 * rOuter), target: new THREE.Vector3(0, 0, 0) };
      },
      readout() {
        const lines = [];
        for (const b of ctx.bodiesNow()) {
          if (b.hidden) continue;
          const p = ctx.planeXYZ(b.sample.x), r = Math.hypot(p.X, p.Y, p.Z);
          if (r >= rP[0] && r <= rP[rP.length - 1]) lines.push(`${b.name}: space falls past at ${ctx.fmt(lookup(vP, r))} c`);
        }
        return lines;
      },
      onOption(id) {
        if (id === 'cut' || id === 'spokes') { rebuild(); return true; }
        return false;
      },
    };
  },
});
