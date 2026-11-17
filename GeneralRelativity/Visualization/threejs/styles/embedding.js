// Embedding style: the equatorial plane as a surface of revolution (Flamm's paraboloid for
// Schwarzschild) coloured by the Kretschmann scalar, reference spheres, scalar field slices, and the
// bodies in 3D. Works for every scene; without a surface it is a plain 3D view.
//
// A compact binary (metric record "binary", Examples.BinaryMerger) gets a time-dependent sheet instead: the
// orbital plane, sunk into wells by the lapse deficit 1 - alpha of the two masses where they are now, and
// rippled by the outgoing gravitational waves. In the orbital plane the leading-order (quadrupole) wave at
// distance r is h(r, phi, t) = (M/r) * 1/2 Re[H(t - r/c) e^{2 i phi}], where H = h+ - i hx is the face-on
// waveform that the scene records (R h / M, the h+ and hx scalars of the tracks). Evaluated at the retarded
// time, it is a two-armed spiral moving out at c, chirping as the orbit shrinks, then the ringdown burst.
(function () {
  const isBinary = (data) => !!(data && data.metric && data.metric.name === 'binary' && data.trackSpecs.length >= 2);

  const VERTEX = /* glsl */`
    attribute vec2 plane;     // (X, Y) in the orbital plane
    attribute float ring;     // index of the ring (radius) of this vertex
    uniform vec2 uX1, uX2;
    uniform float uM1, uM2, uR1, uR2, uMf, uRrem, uRplus;
    uniform int uMerged, uWaves, uFalloff;
    uniform float uDepth, uGain, uRef, uTaperA, uTaperB;
    uniform sampler2D uWave;  // ring k: (Re, Im) of H(t - r_k / c)
    varying vec3 vWorld;
    varying vec2 vPlane;
    varying float vWell, vWave;
    void main() {
      vec2 P = plane;
      float r = length(P);
      float f;
      if (uMerged == 0) {
        float r1 = length(P - uX1), r2 = length(P - uX2);
        // Exterior field of each body; a star's is flat inside its surface.
        f = 2.0 * uM1 / max(r1, max(uR1, 1e-3)) + 2.0 * uM2 / max(r2, max(uR2, 1e-3));
      } else {
        f = 2.0 * uMf / max(r, max(uRrem, 1e-3));
      }
      // Lapse of the (superposed) Kerr-Schild metric, alpha = 1 / sqrt(1 + f); the sheet sinks by 1 - alpha.
      float well = 1.0 - inversesqrt(1.0 + f);
      float y = -uDepth * well;
      float w = 0.0;
      if (uWaves == 1 && r > 0.0) {
        vec2 H = texelFetch(uWave, ivec2(int(ring + 0.5), 0), 0).xy;
        float c2 = (P.x * P.x - P.y * P.y) / (r * r), s2 = 2.0 * P.x * P.y / (r * r);   // cos 2 phi, sin 2 phi
        float q = 0.5 * (H.x * c2 - H.y * s2);                                         // r h / M
        float fall = uFalloff == 1 ? uRef / r : 1.0;
        w = uGain * q * fall * smoothstep(uTaperA, uTaperB, r);
        y += w;
      }
      vWorld = (modelMatrix * vec4(P.x, y, -P.y, 1.0)).xyz;
      vWell = well; vWave = w; vPlane = P;
      gl_Position = projectionMatrix * viewMatrix * vec4(vWorld, 1.0);
    }`;
  const FRAGMENT = /* glsl */`
    uniform float uWaveScale;
    uniform vec2 uX1, uX2;
    uniform float uM1, uM2, uR1, uR2, uRrem, uRplus;
    uniform int uMerged;
    varying vec3 vWorld;
    varying vec2 vPlane;
    varying float vWell, vWave;
    void main() {
      // Inside a horizon (1) or a star (2), decided per pixel so the edges are smooth.
      float inside = 0.0;
      if (uMerged == 0) {
        float r1 = length(vPlane - uX1), r2 = length(vPlane - uX2);
        if (uR1 > 0.0 ? r1 < uR1 : r1 < 2.0 * uM1) inside = uR1 > 0.0 ? 2.0 : 1.0;
        if (uR2 > 0.0 ? r2 < uR2 : r2 < 2.0 * uM2) inside = uR2 > 0.0 ? 2.0 : 1.0;
      } else {
        float r = length(vPlane);
        if (uRrem > 0.0 ? r < uRrem : r < uRplus) inside = uRrem > 0.0 ? 2.0 : 1.0;
      }
      vec3 n = normalize(cross(dFdx(vWorld), dFdy(vWorld)));
      float lit = 0.4 + 0.6 * abs(dot(n, normalize(vec3(0.35, 1.0, 0.45))));
      vec3 base = mix(vec3(0.09, 0.11, 0.17), vec3(0.42, 0.52, 0.85), clamp(vWell / 0.3, 0.0, 1.0));
      float k = clamp(abs(vWave) / max(uWaveScale, 1e-12), 0.0, 1.0);
      vec3 tint = vWave > 0.0 ? vec3(1.0, 0.55, 0.22) : vec3(0.25, 0.62, 1.0);
      vec3 col = mix(base, tint, 0.9 * pow(k, 0.6));
      if (inside > 1.5) col = vec3(0.85, 0.9, 1.0);
      else if (inside > 0.5) col = vec3(0.0);
      gl_FragColor = vec4(col * lit, 1.0);
    }`;
  const LINE_FRAGMENT = /* glsl */`
    varying vec3 vWorld;
    varying vec2 vPlane;
    varying float vWell, vWave;
    void main() { gl_FragColor = vec4(0.75, 0.82, 1.0, 0.16); }`;

  function binarySheet(ctx) {
    const { THREE, data, options } = ctx;
    const group = new THREE.Group();
    const p = data.metric.p || [];
    const bin = { m1: p[0], m2: p[1], R1: p[2] || 0, R2: p[3] || 0, tMerge: p[4] ?? 1e30, Mf: p[5] ?? (p[0] + p[1]), aF: p[6] || 0, Rrem: p[7] || 0 };
    const M = bin.m1 + bin.m2, eta = bin.m1 * bin.m2 / (M * M);
    const rPlus = bin.Rrem > 0 ? 0 : bin.Mf + Math.sqrt(Math.max(bin.Mf * bin.Mf - bin.aF * bin.aF, 0));
    const c = data.header.c || 1;
    const region = options.region * ctx.referenceRadius;
    const [s1, s2] = data.trackSpecs;
    const times = s1.samples.map(ctx.timeOf), times2 = s2.samples.map(ctx.timeOf);

    // Body positions at time t, interpolated in the recorded samples (independent of the core's state).
    const at = (spec, ts, t) => {
      const n = ts.length;
      if (t <= ts[0]) return spec.samples[0].x;
      if (t >= ts[n - 1]) return spec.samples[n - 1].x;
      let lo = 0, hi = n - 1;
      while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (ts[mid] <= t) lo = mid; else hi = mid; }
      const f = (t - ts[lo]) / Math.max(ts[hi] - ts[lo], 1e-30), a = spec.samples[lo].x, b = spec.samples[hi].x;
      return a.map((v, k) => v + f * (b[k] - v));
    };
    let rOrbit = 0;
    for (const sp of [s1, s2]) for (const smp of sp.samples) rOrbit = Math.max(rOrbit, Math.hypot(smp.x[1], smp.x[2]));

    // The face-on waveform H = h+ - i hx (R h / M) from the track scalars; without them, the circular-orbit
    // quadrupole -4 eta (M / d) e^{-2 i Phi} from the positions.
    const iPlus = s1.scalarNames.indexOf('h+'), iCross = s1.scalarNames.indexOf('hx');
    const n = s1.samples.length;
    const Hre = new Float64Array(n), Him = new Float64Array(n);
    s1.samples.forEach((smp, k) => {
      if (iPlus >= 0 && iCross >= 0) { Hre[k] = smp.scalars[iPlus]; Him[k] = -smp.scalars[iCross]; return; }
      const b = at(s2, times2, times[k]), dx = smp.x[1] - b[1], dy = smp.x[2] - b[2], d = Math.hypot(dx, dy);
      const amp = d > 1e-9 ? -4 * eta * M / d : 0, ph = Math.atan2(dy, dx);
      Hre[k] = amp * Math.cos(2 * ph); Him[k] = -amp * Math.sin(2 * ph);
    });
    let qPeak = 0;
    for (let k = 0; k < n; k++) qPeak = Math.max(qPeak, 0.5 * Math.hypot(Hre[k], Him[k]));
    // Before the first sample the binary was already orbiting: continue the first samples' phase rate back in
    // time (otherwise a false wave front would leave the origin at t = 0).
    let k1 = 1;
    while (k1 < n - 1 && times[k1] - times[0] < 8) k1++;
    const arg = Math.atan2(Hre[0] * Him[k1] - Him[0] * Hre[k1], Hre[0] * Hre[k1] + Him[0] * Him[k1]);
    const rate = n > 1 ? arg / Math.max(times[k1] - times[0], 1e-30) : 0;
    const strainAt = (t, out) => {
      if (t <= times[0]) {
        const ph = rate * (t - times[0]), cs = Math.cos(ph), sn = Math.sin(ph);
        out[0] = Hre[0] * cs - Him[0] * sn; out[1] = Hre[0] * sn + Him[0] * cs; return out;
      }
      if (t >= times[n - 1]) { out[0] = Hre[n - 1]; out[1] = Him[n - 1]; return out; }
      let lo = 0, hi = n - 1;
      while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (times[mid] <= t) lo = mid; else hi = mid; }
      const f = (t - times[lo]) / Math.max(times[hi] - times[lo], 1e-30);
      out[0] = Hre[lo] + f * (Hre[hi] - Hre[lo]); out[1] = Him[lo] + f * (Him[hi] - Him[lo]); return out;
    };

    // Polar grid, rings closer together near the bodies (radius ~ k^1.6).
    const Nr = 260, Np = 288;
    const rings = Array.from({ length: Nr }, (_, k) => region * Math.pow(k / (Nr - 1), 1.6));
    const plane = new Float32Array(Nr * Np * 2), ringIdx = new Float32Array(Nr * Np), pos = new Float32Array(Nr * Np * 3);
    for (let k = 0; k < Nr; k++) for (let j = 0; j < Np; j++) {
      const v = k * Np + j, ph = 2 * Math.PI * j / Np, X = rings[k] * Math.cos(ph), Y = rings[k] * Math.sin(ph);
      plane[2 * v] = X; plane[2 * v + 1] = Y; ringIdx[v] = k;
      pos[3 * v] = X; pos[3 * v + 1] = 0; pos[3 * v + 2] = -Y;
    }
    const tri = [], lines = [];
    for (let k = 0; k < Nr - 1; k++) for (let j = 0; j < Np; j++) {
      const a = k * Np + j, b = k * Np + (j + 1) % Np, cc = a + Np, d = b + Np;
      tri.push(a, cc, b, b, cc, d);
      if (j % 12 === 0) lines.push(a, cc);
    }
    for (let k = 0; k < Nr; k += 8) for (let j = 0; j < Np; j++) lines.push(k * Np + j, k * Np + (j + 1) % Np);
    const makeGeometry = (index) => {
      const g = new THREE.BufferGeometry();
      g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
      g.setAttribute('plane', new THREE.BufferAttribute(plane, 2));
      g.setAttribute('ring', new THREE.BufferAttribute(ringIdx, 1));
      g.setIndex(index);
      return g;
    };
    const waveData = new Float32Array(Nr * 4);
    const waveTex = new THREE.DataTexture(waveData, Nr, 1, THREE.RGBAFormat, THREE.FloatType);
    waveTex.minFilter = waveTex.magFilter = THREE.NearestFilter;
    const uniforms = {
      uX1: { value: new THREE.Vector2() }, uX2: { value: new THREE.Vector2() },
      uM1: { value: bin.m1 }, uM2: { value: bin.m2 }, uR1: { value: bin.R1 }, uR2: { value: bin.R2 },
      uMf: { value: bin.Mf }, uRrem: { value: bin.Rrem }, uRplus: { value: rPlus }, uMerged: { value: 0 },
      uWaves: { value: options.waves ? 1 : 0 }, uFalloff: { value: Math.round(options.falloff) },
      uDepth: { value: 0 }, uGain: { value: 0 }, uRef: { value: 0.25 * region },
      uTaperA: { value: 0.5 * rOrbit }, uTaperB: { value: 1.6 * rOrbit }, uWave: { value: waveTex }, uWaveScale: { value: 1 },
    };
    const sheet = new THREE.Mesh(makeGeometry(tri), new THREE.ShaderMaterial({ vertexShader: VERTEX, fragmentShader: FRAGMENT, uniforms, side: THREE.DoubleSide }));
    const wire = new THREE.LineSegments(makeGeometry(lines), new THREE.ShaderMaterial({ vertexShader: VERTEX, fragmentShader: LINE_FRAGMENT, uniforms, transparent: true, depthWrite: false }));
    sheet.frustumCulled = wire.frustumCulled = false;   // displaced on the GPU
    sheet.visible = !!options.surface;
    wire.visible = !!options.surface && !!options.wire;
    group.add(sheet, wire);

    // Display scales: wells of depth `depth` x 6 M at the centre; the strongest wave (the merger) rises
    // `waveGain` x 6% of the region.
    const setScales = () => {
      uniforms.uDepth.value = options.depth * 6 * M;
      uniforms.uGain.value = qPeak > 0 ? options.waveGain * 0.06 * region / qPeak : 0;
      uniforms.uWaveScale.value = uniforms.uGain.value * qPeak * 0.7;
      uniforms.uWaves.value = options.waves ? 1 : 0;
      uniforms.uFalloff.value = Math.round(options.falloff);
    };
    setScales();
    const diverging = (x) => {
      const t = Math.min(1, Math.max(0, x)), u = Math.abs(2 * t - 1);
      const lo = new THREE.Color(0.25, 0.62, 1.0), hi = new THREE.Color(1.0, 0.55, 0.22), mid = new THREE.Color(0.09, 0.11, 0.17);
      return mid.clone().lerp(t < 0.5 ? lo : hi, Math.pow(u, 0.6) * 0.9);
    };
    ctx.drawScale(-qPeak, qPeak, 'r h / M in the orbital plane', diverging);
    ctx.setNotes([
      'Sheet: the orbital plane. Depth: the clocks\' lag 1 - alpha of the superposed Kerr-Schild metric of the two masses where they are now (the remnant after the merger); black: inside a horizon.',
      'Ripples: the outgoing gravitational wave in the orbital plane, h = (M/r) 1/2 Re[(h+ - i hx)(t - r/c) e^{2 i phi}], from the face-on waveform the simulation recorded, at the retarded time: a two-armed spiral moving out at c. Orange: stretched along the line of sight to the source, blue: squeezed.',
      `Exaggerated for display: ${options.falloff ? 'true 1/r fall-off' : 'shown as r h (no 1/r fall-off)'}, the merger wave rises ${ctx.fmt(options.waveGain * 0.06 * region)} M; the near zone (r < ${ctx.fmt(1.6 * rOrbit)}) is faded out, where the wave formula does not apply.`,
    ].join('\n'));

    // State at the clock: body positions, merged or not, and H at the retarded time of every ring.
    let stateT = NaN;
    const P1 = [0, 0], P2 = [0, 0], tmp = [0, 0];
    let merged = false;
    const stateAt = (t) => {
      if (t === stateT) return;
      stateT = t;
      const a = at(s1, times, t), b = at(s2, times2, t);
      P1[0] = a[1]; P1[1] = a[2]; P2[0] = b[1]; P2[1] = b[2];
      merged = c * t >= bin.tMerge;
      uniforms.uX1.value.set(P1[0], P1[1]);
      uniforms.uX2.value.set(P2[0], P2[1]);
      uniforms.uMerged.value = merged ? 1 : 0;
      for (let k = 0; k < Nr; k++) {
        strainAt(t - rings[k] / c, tmp);
        waveData[4 * k] = tmp[0]; waveData[4 * k + 1] = tmp[1];
      }
      waveTex.needsUpdate = true;
    };
    // The same height as the vertex shader, for the bodies and their paths.
    const heightAt = (X, Y, t) => {
      stateAt(t);
      const r = Math.hypot(X, Y);
      let f;
      if (!merged) {
        f = 2 * bin.m1 / Math.max(Math.hypot(X - P1[0], Y - P1[1]), bin.R1, 1e-3) + 2 * bin.m2 / Math.max(Math.hypot(X - P2[0], Y - P2[1]), bin.R2, 1e-3);
      } else f = 2 * bin.Mf / Math.max(r, bin.Rrem, 1e-3);
      let y = -uniforms.uDepth.value * (1 - 1 / Math.sqrt(1 + f));
      if (options.waves && r > 0) {
        strainAt(t - r / c, tmp);
        const q = 0.5 * (tmp[0] * (X * X - Y * Y) - tmp[1] * 2 * X * Y) / (r * r);
        const fall = options.falloff ? uniforms.uRef.value / r : 1;
        const e = Math.min(1, Math.max(0, (r - uniforms.uTaperA.value) / (uniforms.uTaperB.value - uniforms.uTaperA.value)));
        y += uniforms.uGain.value * q * fall * e * e * (3 - 2 * e);
      }
      return y;
    };

    return {
      group,
      dynamicPlacement: true,
      bodyRadius: () => Math.max(0.35, 0.004 * region),
      place(x, sample, t) {
        const time = t === undefined ? ctx.clock : t;
        return ctx.toWorld(x[1], x[2], heightAt(x[1], x[2], time) + 0.15);
      },
      update(clock) { stateAt(clock); },
      camera() {
        return { position: ctx.toWorld(0, -1.0 * region, 0.62 * region), target: new THREE.Vector3(0, -0.05 * region, 0) };
      },
      readout() {
        const t = ctx.clock, lines = [];
        lines.push(merged ? `merged at ct = ${ctx.fmt(bin.tMerge)}; its signal has travelled to r = ${ctx.fmt(Math.max(0, c * t - bin.tMerge))}` : `inspiral: separation ${ctx.fmt(Math.hypot(P1[0] - P2[0], P1[1] - P2[1]))}`);
        strainAt(t - 0.5 * region / c, tmp);
        lines.push(`at r = ${ctx.fmt(0.5 * region)}: r h / M up to ${ctx.fmt(0.5 * Math.hypot(tmp[0], tmp[1]))} (peak ${ctx.fmt(qPeak)})`);
        return lines;
      },
      onOption(id, v) {
        if (id === 'surface' || id === 'wire') { sheet.visible = !!options.surface; wire.visible = !!options.surface && !!options.wire; return true; }
        if (id === 'depth' || id === 'waveGain' || id === 'waves' || id === 'falloff') { setScales(); stateT = NaN; return false; }
        return false;
      },
      dispose() { waveTex.dispose(); },
    };
  }

(window.GRStyles = window.GRStyles || []).push({
  id: 'embedding',
  label: 'Embedding / 3D',
  order: 1,
  description: 'The equatorial plane embedded in Euclidean space (depth shows the stretching of radial distance), coloured by curvature; spheres and field slices in 3D. A merging binary: the orbital plane with the wells of the two masses and the outgoing gravitational waves.',
  options: [
    { id: 'surface', label: 'surface', type: 'checkbox', value: 1 },
    { id: 'wire', label: 'wireframe', type: 'checkbox', value: 1 },
    { id: 'fields', label: 'field slices', type: 'checkbox', value: 1, show: (data) => !isBinary(data) },
    { id: 'waves', label: 'gravitational waves', type: 'checkbox', value: 1, show: isBinary },
    { id: 'waveGain', label: 'wave height', type: 'range', min: 0, max: 6, step: 0.1, value: 1.5, show: isBinary, help: 'Exaggeration: the merger wave rises this times 6% of the region' },
    { id: 'falloff', label: 'amplitude', type: 'select', value: 0, show: isBinary,
      choices: [{ value: 0, label: 'r h (no fall-off)' }, { value: 1, label: 'h (true 1/r)' }] },
    { id: 'region', label: 'region (M)', type: 'range', min: 20, max: 400, step: 5, value: 100, show: isBinary, help: 'Radius of the sheet; the waves need room: their wavelength is 160 M at the start and 18 M at the merger' },
    { id: 'depth', label: 'well depth', type: 'range', min: 0, max: 6, step: 0.1, value: 2.5, show: isBinary, help: 'The sheet sinks by this times 6 M times 1 - alpha' },
  ],
  create(ctx) {
    if (isBinary(ctx.data)) return binarySheet(ctx);
    const { THREE, data, options } = ctx;
    const group = new THREE.Group();
    const rows = data.rows;
    let surfaceAt = null;   // r -> {R, z}
    const surfaceObjects = [];
    const fieldObjects = [];

    if (rows.length) {
      const nPhi = rows[0].nPhi;
      const [lo, hi] = ctx.logRange(rows.map(r => r.kretschmann));
      ctx.drawScale(lo, hi, 'log10 Kretschmann');
      const rs = rows.map(r => r.r), Rs = rows.map(r => r.R ?? r.r), zs = rows.map(r => r.z);
      const interp = (arr, r) => {
        if (r <= rs[0]) return arr[0];
        if (r >= rs[rs.length - 1]) return arr[arr.length - 1];
        let i = 1;
        while (rs[i] < r) i++;
        const f = (r - rs[i - 1]) / (rs[i] - rs[i - 1]);
        return arr[i - 1] + f * (arr[i] - arr[i - 1]);
      };
      surfaceAt = (r) => ({ R: r >= rs[rs.length - 1] ? r * Rs[Rs.length - 1] / rs[rs.length - 1] : interp(Rs, r), z: interp(zs, r) });
      const positions = [], colors = [], indices = [];
      rows.forEach((row) => {
        const t = hi > lo ? (Math.log10(Math.max(row.kretschmann, 1e-300)) - lo) / (hi - lo) : 0.5;
        const col = ctx.colormap(t);
        const R = row.R ?? row.r;
        for (let j = 0; j < nPhi; j++) {
          const ph = 2 * Math.PI * j / nPhi;
          positions.push(R * Math.cos(ph), row.z, -R * Math.sin(ph));
          colors.push(col.r, col.g, col.b);
        }
      });
      for (let i = 0; i < rows.length - 1; i++) {
        for (let j = 0; j < nPhi; j++) {
          const a = i * nPhi + j, b = i * nPhi + (j + 1) % nPhi, c = a + nPhi, d = b + nPhi;
          indices.push(a, c, b, b, c, d);
        }
      }
      const geom = new THREE.BufferGeometry();
      geom.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
      geom.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
      geom.setIndex(indices);
      geom.computeVertexNormals();
      const mesh = new THREE.Mesh(geom, new THREE.MeshStandardMaterial({ vertexColors: true, side: THREE.DoubleSide, roughness: 0.8, metalness: 0.05 }));
      const wire = new THREE.LineSegments(new THREE.WireframeGeometry(geom), new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.08 }));
      mesh.visible = !!options.surface;
      wire.visible = !!options.surface && !!options.wire;
      group.add(mesh, wire);
      surfaceObjects.push(mesh, wire);
    }

    for (const field of data.fields) {
      const nu = field.u.length, nv = field.v.length;
      const [lo, hi] = ctx.logRange(field.values.flat());
      if (!rows.length) ctx.drawScale(lo, hi, 'log10 ' + field.name);
      const positions = [], colors = [], indices = [];
      for (let i = 0; i < nu; i++) {
        for (let j = 0; j < nv; j++) {
          const u = field.u[i], v = field.v[j];
          const p = field.plane === 'xy' ? [u, v, field.offset] : field.plane === 'xz' ? [u, field.offset, v] : [field.offset, u, v];
          positions.push(p[0], p[2], -p[1]);
          const val = field.values[i][j];
          const col = ctx.colormap(hi > lo ? (Math.log10(Math.max(val, 1e-300)) - lo) / (hi - lo) : 0.5);
          colors.push(col.r, col.g, col.b);
        }
      }
      for (let i = 0; i < nu - 1; i++) {
        for (let j = 0; j < nv - 1; j++) {
          const a = i * nv + j, b = a + 1, c = a + nv, d = c + 1;
          indices.push(a, c, b, b, c, d);
        }
      }
      const geom = new THREE.BufferGeometry();
      geom.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
      geom.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
      geom.setIndex(indices);
      const mesh = new THREE.Mesh(geom, new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide, transparent: true, opacity: 0.85 }));
      mesh.visible = !!options.fields;
      group.add(mesh);
      fieldObjects.push(mesh);
    }

    for (const s of data.spheres) {
      const mesh = new THREE.Mesh(new THREE.SphereGeometry(s.radius, 48, 32),
        new THREE.MeshStandardMaterial({ color: new THREE.Color(s.color || '#888888'), transparent: true, opacity: 0.55, roughness: 0.6 }));
      const c = s.centre || [0, 0, 0];
      mesh.position.copy(ctx.toWorld(c[0], c[1], c[2]));
      group.add(mesh);
    }

    const lift = 0.02 * ctx.referenceRadius;
    return {
      group,
      place(x) {
        if (data.header.chart === 'spherical' && surfaceAt && Math.abs(x[2] - Math.PI / 2) < 1e-6) {
          const s = surfaceAt(x[1]);
          return new THREE.Vector3(s.R * Math.cos(x[3]), s.z + lift, -s.R * Math.sin(x[3]));
        }
        const p = ctx.planeXYZ(x);
        return ctx.toWorld(p.X, p.Y, p.Z);
      },
      onOption(id, v) {
        if (id === 'surface' || id === 'wire') {
          if (surfaceObjects[0]) surfaceObjects[0].visible = !!options.surface;
          if (surfaceObjects[1]) surfaceObjects[1].visible = !!options.surface && !!options.wire;
          return true;
        }
        if (id === 'fields') { fieldObjects.forEach(m => { m.visible = !!v; }); return true; }
        return false;
      },
    };
  },
});
})();
