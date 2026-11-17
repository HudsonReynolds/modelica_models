// Lensing style: a telescope-like image of background galaxies seen through a foreground lens, in the
// thin-lens approximation. Per pixel at angle theta the shader solves the lens equation
// beta = theta - alpha(theta) and shades the unlensed sky at beta; lensing conserves surface brightness,
// so arcs, rings and multiple images appear by themselves. Lens models and cosmological distances are
// the formulas of GeneralRelativity.Lensing (tested there): point mass, singular isothermal ellipse with
// external shear, flat Lambda-CDM angular-diameter distances. Angles are in arcseconds.
(function () {
  const C = 299792458.0, G = 6.67430e-11, GM_SUN = 1.32712440018e20, MPC = 3.0856775814913673e22;
  const ARCSEC = Math.PI / 648000;

  // Flat-universe angular-diameter distance (Simpson's rule), as Lensing.angularDiameterDistance.
  function angularDiameterDistance(z1, z2, H0, Om) {
    const DH = C / (H0 * 1000 / MPC), n = 400;
    const chi = (z) => {
      const h = z / n; let s = 0;
      for (let i = 0; i <= n; i++) {
        const zi = i * h, w = (i === 0 || i === n) ? 1 : (i % 2 ? 4 : 2);
        s += w / Math.sqrt(Om * Math.pow(1 + zi, 3) + 1 - Om);
      }
      return DH * s * h / 3;
    };
    return (chi(z2) - chi(z1)) / (1 + z2);
  }
  // The lens in physical terms, from the sliders.
  function lensPhysics(o, H0, Om) {
    const zl = o.zl, zs = Math.max(o.zs, o.zl + 0.01);
    const Dl = angularDiameterDistance(0, zl, H0, Om), Ds = angularDiameterDistance(0, zs, H0, Om), Dls = angularDiameterDistance(zl, zs, H0, Om);
    const GM = Math.pow(10, o.mass) * GM_SUN;
    const aperture = o.aperture * 1e-3 * MPC;
    const sigma = Math.sqrt(GM / (Math.PI * aperture));
    const point = o.model === 1;
    const thetaE = point ? Math.sqrt(4 * GM / (C * C) * Dls / (Dl * Ds)) : 4 * Math.PI * Math.pow(sigma / C, 2) * Dls / Ds;
    const SigmaCr = C * C * Ds / (4 * Math.PI * G * Dl * Dls);
    return { zl, zs, Dl, Ds, Dls, sigma, thetaE: thetaE / ARCSEC, SigmaCr, massInRing: Math.PI * Math.pow(Dl * thetaE, 2) * SigmaCr * G / GM_SUN, point };
  }
  // Scaled deflection (arcsec) at image position (x, y) for the current lens, as Lensing.deflection*.
  function deflection(x, y, L) {
    let dx = x - L.lx, dy = y - L.ly;
    let ax = 0, ay = 0;
    if (L.point) {
      const r2 = dx * dx + dy * dy;
      if (r2 > 0) { ax = L.b * L.b * dx / r2; ay = L.b * L.b * dy / r2; }
    } else {
      const c = Math.cos(L.pa), s = Math.sin(L.pa);
      const u = c * dx + s * dy, v = -s * dx + c * dy;
      const q = L.q, sq = Math.sqrt(Math.max(1 - q * q, 0)), psi = Math.sqrt(q * q * u * u + v * v);
      let au = 0, av = 0;
      if (psi > 0) {
        if (sq < 1e-6) { au = L.b * u / psi; av = L.b * v / psi; }
        else { au = L.b * q / sq * Math.atan(sq * u / psi); av = L.b * q / sq * Math.atanh(Math.min(sq * v / psi, 0.999999999)); }
      }
      ax = c * au - s * av; ay = s * au + c * av;
    }
    const g1 = L.gamma * Math.cos(2 * L.gpa), g2 = L.gamma * Math.sin(2 * L.gpa);
    ax += g1 * dx + g2 * dy; ay += g2 * dx - g1 * dy;
    return [ax, ay];
  }
  // Seeded background galaxies, roughly in deep-field proportions: clumpy star-forming irregulars, spirals
  // and ellipticals, mostly small. One galaxy is placed within 0.6" of the field centre so that a lens there
  // shows arcs from the start; it is otherwise drawn like the rest.
  function makeSources(n, seed, spread) {
    let state = (seed * 2654435761) >>> 0;
    const rand = () => { state = (state * 1664525 + 1013904223) >>> 0; return state / 4294967296; };
    const out = [];
    for (let i = 0; i < n; i++) {
      const u = rand();
      const type = i === 0 ? 1 : (u < 0.42 ? 2 : u < 0.8 ? 1 : 0);
      const size = (i === 0 ? 0.16 : 0.05) * Math.exp(0.55 * Math.sqrt(-2 * Math.log(rand() + 1e-9)) * Math.cos(2 * Math.PI * rand()));
      const near = i === 0;
      const a = 2 * Math.PI * rand(), rr = 0.6 * Math.sqrt(rand());
      out.push({
        x: near ? rr * Math.cos(a) : (rand() * 2 - 1) * spread, y: near ? rr * Math.sin(a) : (rand() * 2 - 1) * spread,
        size: Math.min(Math.max(size, 0.02), 0.6),
        q: type === 0 ? 0.55 + 0.45 * rand() : 0.25 + 0.75 * rand(), pa: rand() * Math.PI,
        bright: (type === 2 ? 0.5 : 0.8) + 0.8 * rand(), type, hue: rand(),
        pitch: 0.3 + 0.3 * rand(), arms: rand() < 0.75 ? 2 : 3, bulge: 0.3 + 0.7 * rand(), phase: 2 * Math.PI * rand(),
      });
    }
    return out;
  }
  // Foreground stars (in our own galaxy, so not lensed), with diffraction spikes.
  function makeStars(seed, spread) {
    let state = (seed * 40503 + 17) >>> 0;
    const rand = () => { state = (state * 1664525 + 1013904223) >>> 0; return state / 4294967296; };
    return Array.from({ length: 3 }, () => ({ x: (rand() * 2 - 1) * spread, y: (rand() * 2 - 1) * spread, bright: 0.5 + 2.5 * rand() * rand() }));
  }

  const VERTEX = 'varying vec2 vUv; void main() { vUv = uv; gl_Position = vec4(position.xy, 0.0, 1.0); }';
  // Source parameters live in a float texture (3 texels per galaxy), not uniforms, so the count is not limited
  // by the GPU's uniform budget. A Hubble deep field has ~0.25 galaxies per square arcsec, ~150 in the 24" sky.
  const MAXSRC = 400;
  // Spatial index: the source plane is cut into GRID x GRID cells, and each cell lists the galaxies whose
  // footprint (radius 7 scale lengths, where the profiles are cut off) overlaps it. A pixel then visits only
  // the few galaxies in the cell its ray lands in, instead of all of them.
  const GRID = 64, IDXW = 256;
  const FRAGMENT = /* glsl */`
    varying vec2 vUv;
    uniform vec2 uCentre; uniform float uFov, uAspect;
    uniform int uLensing, uPoint, uCrit, uMap, uSuper;
    uniform float uB, uQ, uPA, uGamma, uGPA, uLensLight, uExposure, uRing, uPixel;
    uniform vec2 uLens;
    uniform int uCount;
    // Row i: texel 0 = x, y, size, pa; 1 = q, brightness, type (0 elliptical, 1 spiral, 2 irregular), hue;
    // 2 = pitch (tan), arms, bulge fraction, arm phase.
    uniform highp sampler2D uSrc;
    uniform highp sampler2D uCells;  // GRID x GRID: first list entry, galaxy count
    uniform highp sampler2D uIndex;  // the cells' galaxy lists, IDXW entries per row
    uniform float uGridSpan;         // the cells cover [-uGridSpan, uGridSpan]^2 (arcsec)
    uniform int uStars;
    uniform vec4 uStar[4];           // x, y, brightness, unused
    float hash2(vec2 p) { return fract(sin(dot(p, vec2(127.1, 311.7))) * 43758.5453); }
    float vnoise(vec2 p) {
      vec2 i = floor(p), f = fract(p);
      vec2 u = f * f * (3.0 - 2.0 * f);
      return mix(mix(hash2(i), hash2(i + vec2(1.0, 0.0)), u.x), mix(hash2(i + vec2(0.0, 1.0)), hash2(i + vec2(1.0, 1.0)), u.x), u.y);
    }
    float fbm(vec2 p) { return 0.5 * vnoise(p) + 0.3 * vnoise(2.03 * p + 7.1) + 0.2 * vnoise(4.07 * p + 3.3); }
    vec2 deflect(vec2 th) {
      vec2 d = th - uLens;
      vec2 a = vec2(0.0);
      if (uPoint == 1) {
        float r2 = dot(d, d);
        a = r2 > 0.0 ? uB * uB * d / r2 : vec2(0.0);
      } else {
        float c = cos(uPA), s = sin(uPA);
        vec2 w = vec2(c * d.x + s * d.y, -s * d.x + c * d.y);
        float sq = sqrt(max(1.0 - uQ * uQ, 0.0));
        float psi = sqrt(uQ * uQ * w.x * w.x + w.y * w.y);
        vec2 aw = vec2(0.0);
        if (psi > 0.0) {
          if (sq < 1e-4) aw = uB * w / psi;
          else {
            float t = clamp(sq * w.y / psi, -0.9999999, 0.9999999);
            aw = uB * uQ / sq * vec2(atan(sq * w.x / psi), 0.5 * log((1.0 + t) / (1.0 - t)));
          }
        }
        a = vec2(c * aw.x - s * aw.y, s * aw.x + c * aw.y);
      }
      float g1 = uGamma * cos(2.0 * uGPA), g2 = uGamma * sin(2.0 * uGPA);
      return a + vec2(g1 * d.x + g2 * d.y, g2 * d.x - g1 * d.y);
    }
    // One galaxy at source-plane position beta. Coordinates are deprojected into the galaxy's own plane
    // (rotate by the position angle, stretch the minor axis by 1/q) in units of its scale length.
    vec3 galaxy(vec2 beta, vec4 A, vec4 B, vec4 Cc) {
      vec2 d = beta - A.xy;
      if (dot(d, d) > 49.0 * A.z * A.z) return vec3(0.0);
      float c = cos(A.w), s = sin(A.w);
      vec2 w = vec2(c * d.x + s * d.y, (-s * d.x + c * d.y) / B.x) / A.z;
      float rr = length(w);
      if (rr > 7.0) return vec3(0.0);
      vec2 seed = vec2(B.w * 91.7, B.w * 37.3);
      if (B.z < 0.5) {
        // Elliptical: de Vaucouleurs profile (effective radius 1), old red-orange stars, no structure.
        float I = exp(-7.67 * (pow(rr + 1e-3, 0.25) - 1.0));
        return B.y * 0.004 * I * mix(vec3(1.0, 0.68, 0.42), vec3(1.0, 0.8, 0.58), B.w);
      }
      if (B.z > 1.5) {
        // Irregular, star-forming: an exponential envelope broken into bright blue clumps.
        float clumps = pow(max(fbm(w * 2.2 + seed), 0.0), 3.0) * 4.0;
        float I = exp(-1.2 * rr) * (0.25 + clumps);
        return B.y * 0.22 * I * mix(vec3(0.55, 0.7, 1.0), vec3(0.75, 0.75, 1.0), B.w);
      }
      // Spiral: Sersic bulge + exponential disk; logarithmic arms (pitch tan = Cc.x) that start outside the
      // bulge; star-forming knots along the arms; dust lanes on the arms' inner edges.
      float ph = atan(w.y, w.x);
      float arg = Cc.y * (ph - log(rr + 1e-3) / Cc.x - Cc.w);
      // max(0, .): rounding can push 0.5 + 0.5 cos slightly negative, and pow of a negative base is NaN.
      float arm = pow(max(0.5 + 0.5 * cos(arg), 0.0), 2.5);
      float dust = pow(max(0.5 + 0.5 * cos(arg + 0.8), 0.0), 8.0);
      float ramp = smoothstep(0.35, 0.9, rr) * (1.0 - smoothstep(2.6, 4.5, rr));
      float knots = pow(max(fbm(w * 5.0 + seed), 0.0), 2.5) * 2.0;
      float disk = exp(-rr) * (0.7 + 1.1 * arm * ramp * (0.7 + knots)) * (1.0 - 0.45 * dust * ramp);
      float bulge = Cc.z * 6.0 * exp(-3.67 * (pow(rr / 0.25 + 1e-3, 0.5) - 1.0)) * 0.02;
      vec3 diskColour = mix(vec3(1.0, 0.9, 0.75), vec3(0.6, 0.72, 1.0), smoothstep(0.3, 1.8, rr));
      vec3 knotColour = mix(vec3(0.65, 0.78, 1.0), vec3(1.0, 0.65, 0.8), step(0.8, fbm(w * 9.0 + seed.yx)));
      return B.y * 0.25 * (disk * diskColour + 0.5 * arm * ramp * knots * exp(-rr) * knotColour + bulge * vec3(1.0, 0.82, 0.55));
    }
    vec3 sky(vec2 beta) {
      vec3 col = vec3(0.0);
      vec2 g = (beta + uGridSpan) / (2.0 * uGridSpan) * ${GRID}.0;
      // No galaxy reaches outside the grid.
      if (uCount == 0 || any(lessThan(g, vec2(0.0))) || any(greaterThanEqual(g, vec2(${GRID}.0)))) return col;
      vec4 cell = texelFetch(uCells, ivec2(g), 0);
      int start = int(cell.x + 0.5), count = int(cell.y + 0.5);
      for (int k = 0; k < ${MAXSRC}; k++) {
        if (k >= count) break;
        int j = start + k;
        int i = int(texelFetch(uIndex, ivec2(j - ${IDXW} * (j / ${IDXW}), j / ${IDXW}), 0).x + 0.5);
        col += galaxy(beta, texelFetch(uSrc, ivec2(0, i), 0), texelFetch(uSrc, ivec2(1, i), 0), texelFetch(uSrc, ivec2(2, i), 0));
      }
      return col;
    }
    vec3 shade(vec2 th) {
      vec2 beta = uLensing == 1 ? th - deflect(th) : th;
      vec3 col = sky(beta);
      if (uLensing == 1 && uStars == 1) {
        // Foreground stars: a compact core and four diffraction spikes (in front of the lens, not lensed).
        for (int i = 0; i < 4; i++) {
          if (i >= 3) break;
          vec2 d = th - uStar[i].xy;
          float r = length(d);
          float core = exp(-r * r / (2.0 * 0.012 * 0.012)) + 0.15 * exp(-r / 0.05);
          float spikes = exp(-abs(d.x) / 0.004) * exp(-abs(d.y) / 0.35) + exp(-abs(d.y) / 0.004) * exp(-abs(d.x) / 0.35);
          col += uStar[i].z * vec3(1.0, 0.95, 0.9) * (2.0 * core + 0.25 * spikes);
        }
      }
      if (uLensing == 1 && uLensLight > 0.0) {
        // The lens galaxy's own light (de Vaucouleurs profile, effective radius ~ the Einstein radius).
        vec2 d = th - uLens;
        float c = cos(uPA), s = sin(uPA);
        vec2 w = vec2(c * d.x + s * d.y, (-s * d.x + c * d.y) / (uPoint == 1 ? 1.0 : uQ));
        float r = length(w) / max(uB, 0.3);
        col += uLensLight * vec3(1.0, 0.78, 0.5) * 2.5 * exp(-7.67 * (pow(r + 1e-4, 0.25) - 1.0)) * 0.02;
      }
      return col;
    }
    void main() {
      vec2 th0 = uCentre + (vUv * 2.0 - 1.0) * vec2(uFov * uAspect, uFov);
      vec3 col = vec3(0.0);
      // 2 x 2 supersampling (one pixel = uPixel arcsec) against aliasing in thin arcs.
      if (uSuper == 1) {
        for (int k = 0; k < 4; k++) {
          vec2 o = (vec2(float(k - 2 * (k / 2)), float(k / 2)) - 0.5) * 0.5 * uPixel;
          col += 0.25 * shade(th0 + o);
        }
      } else col = shade(th0);
      // Sky background with faint fixed-pattern noise, like a long exposure.
      float n = fract(sin(dot(floor(gl_FragCoord.xy), vec2(12.9898, 78.233))) * 43758.5453);
      col += vec3(0.004, 0.005, 0.008) * (0.6 + 0.8 * n);
      col = vec3(1.0) - exp(-col * uExposure);
      col = pow(col, vec3(1.0 / 2.2));
      if (uLensing == 1) {
        // Inverse magnification det A = (1 - da_x/dx)(1 - da_y/dy) - (da_x/dy)(da_y/dx), by central differences.
        float h = 0.5 * uPixel;
        vec2 ax1 = deflect(th0 + vec2(h, 0.0)), ax0 = deflect(th0 - vec2(h, 0.0));
        vec2 ay1 = deflect(th0 + vec2(0.0, h)), ay0 = deflect(th0 - vec2(0.0, h));
        vec2 dax = (ax1 - ax0) / (2.0 * h), day = (ay1 - ay0) / (2.0 * h);
        float det = (1.0 - dax.x) * (1.0 - day.y) - day.x * dax.y;
        float kappa = 0.5 * (dax.x + day.y);
        if (uMap == 1) col = mix(col, vec3(0.2, 0.6, 1.0) * clamp(kappa, 0.0, 2.0) * 0.6 + col * 0.4, 0.55);
        if (uMap == 2) {
          float lm = clamp(log(1.0 / max(abs(det), 1e-4)) / log(100.0), 0.0, 1.0);
          col = mix(col, mix(vec3(0.1, 0.1, 0.3), vec3(1.0, 0.4, 0.2), lm), 0.5);
        }
        // The singular centre flips the sign of det A without being a critical curve: mask it.
        if (uCrit == 1 && length(th0 - uLens) > 0.03 * max(uB, 0.05)) {
          float w = fwidth(det);
          col = mix(col, vec3(1.0, 0.3, 0.3), (1.0 - smoothstep(0.0, 1.5 * w, abs(det))) * 0.9);
        }
        if (uRing > 0.0) {
          float rr = length(th0 - uLens);
          float w = fwidth(rr);
          col = mix(col, vec3(0.4, 0.9, 1.0), (1.0 - smoothstep(0.0, 1.2 * w, abs(rr - uRing))) * 0.6);
        }
      }
      gl_FragColor = vec4(col, 1.0);
    }`;

  window.GRStyles = window.GRStyles || [];
  window.GRStyles.push({
    id: 'lensing',
    label: 'Gravitational lensing',
    order: 7,
    description: 'Background galaxies seen through a foreground galaxy or point mass, as in a deep telescope image. Change the lens mass, shape and redshifts and watch arcs, Einstein rings and multiple images form. Drag to move the lens; scroll to zoom.',
    options: [
      { id: 'mass', label: 'mass log10 M/M☉', type: 'range', min: 9, max: 13.5, step: 0.02, value: 11.7, help: 'Point: total mass. Galaxy: projected mass inside the aperture.' },
      { id: 'model', label: 'lens', type: 'select', choices: [{ value: 0, label: 'galaxy (isothermal ellipse)' }, { value: 1, label: 'point mass' }], value: 0 },
      { id: 'q', label: 'axis ratio', type: 'range', min: 0.3, max: 1, step: 0.01, value: 0.75 },
      { id: 'pa', label: 'position angle', type: 'range', min: 0, max: 180, step: 1, value: 30 },
      { id: 'gamma', label: 'external shear', type: 'range', min: 0, max: 0.3, step: 0.005, value: 0.05 },
      { id: 'gpa', label: 'shear angle', type: 'range', min: 0, max: 180, step: 1, value: 110 },
      { id: 'zl', label: 'lens redshift', type: 'range', min: 0.05, max: 3, step: 0.01, value: 0.5 },
      { id: 'zs', label: 'source redshift', type: 'range', min: 0.1, max: 6, step: 0.01, value: 2 },
      { id: 'aperture', label: 'aperture (kpc)', type: 'range', min: 1, max: 50, step: 1, value: 10, help: 'Radius that the galaxy mass refers to' },
      { id: 'sources', label: 'galaxies', type: 'range', min: 0, max: MAXSRC, step: 1, value: 150 },
      { id: 'stars', label: 'foreground stars', type: 'checkbox', value: 1 },
      { id: 'seed', label: 'sky seed', type: 'range', min: 1, max: 200, step: 1, value: 7 },
      { id: 'fov', label: 'field (half, ")', type: 'range', min: 0.5, max: 40, step: 0.1, value: 4 },
      { id: 'lensLight', label: 'lens light', type: 'range', min: 0, max: 1.5, step: 0.05, value: 0.6 },
      { id: 'crit', label: 'critical curves', type: 'checkbox', value: 1 },
      { id: 'ring', label: 'Einstein radius', type: 'checkbox', value: 0 },
      { id: 'map', label: 'overlay', type: 'select', choices: [{ value: 0, label: 'none' }, { value: 1, label: 'convergence' }, { value: 2, label: 'magnification' }], value: 0 },
      { id: 'inset', label: 'source plane', type: 'checkbox', value: 1, help: 'Inset: the unlensed sky with the caustics' },
      { id: 'supersample', label: 'anti-alias', type: 'checkbox', value: 1 },
      { id: 'exposure', label: 'exposure', type: 'range', min: 0.2, max: 6, step: 0.1, value: 1.6 },
    ],
    available(data) {
      return data.records && data.records.lens ? null : 'needs a lens scene (Examples.GalaxyLensing or Visualization.writeLens)';
    },
    create(ctx) {
      const { THREE, data, options, renderer, controls } = ctx;
      const scene = data.records.lens[0];
      const H0 = scene.H0 || 70, Om = scene.OmegaM || 0.3;
      // First time: take the lens from the Modelica scene.
      if (!options._fromScene) {
        Object.assign(options, {
          mass: Math.log10(scene.massSolar), model: scene.model === 'point' ? 1 : 0, q: scene.axisRatio, pa: scene.positionAngleDeg,
          gamma: scene.shear, gpa: scene.shearAngleDeg, zl: scene.zLens, zs: scene.zSource, aperture: scene.apertureKpc, _fromScene: true,
        });
      }
      controls.enabled = false;
      const srcData = new Float32Array(3 * MAXSRC * 4);
      const srcTexture = new THREE.DataTexture(srcData, 3, MAXSRC, THREE.RGBAFormat, THREE.FloatType);
      srcTexture.minFilter = srcTexture.magFilter = THREE.NearestFilter;
      const cellData = new Float32Array(GRID * GRID * 4);
      const cellTexture = new THREE.DataTexture(cellData, GRID, GRID, THREE.RGBAFormat, THREE.FloatType);
      cellTexture.minFilter = cellTexture.magFilter = THREE.NearestFilter;
      let indexTexture = null;
      const uniforms = {
        uCentre: { value: new THREE.Vector2() }, uFov: { value: 4 }, uAspect: { value: 1 },
        uLensing: { value: 1 }, uPoint: { value: 0 }, uCrit: { value: 1 }, uMap: { value: 0 }, uSuper: { value: 1 },
        uB: { value: 1 }, uQ: { value: 0.75 }, uPA: { value: 0 }, uGamma: { value: 0 }, uGPA: { value: 0 },
        uLensLight: { value: 0.6 }, uExposure: { value: 1.6 }, uRing: { value: 0 }, uPixel: { value: 0.01 },
        uLens: { value: new THREE.Vector2() }, uCount: { value: 0 },
        uSrc: { value: srcTexture }, uCells: { value: cellTexture }, uIndex: { value: null }, uGridSpan: { value: 16 },
        uStars: { value: 1 }, uStar: { value: Array.from({ length: 4 }, () => new THREE.Vector4()) },
      };
      const material = new THREE.ShaderMaterial({ vertexShader: VERTEX, fragmentShader: FRAGMENT, uniforms, depthTest: false, depthWrite: false });
      const quadScene = new THREE.Scene();
      quadScene.add(new THREE.Mesh(new THREE.PlaneGeometry(2, 2), material));
      const quadCamera = new THREE.OrthographicCamera(-1, 1, 1, -1, 0, 1);
      // Caustic overlay for the source-plane inset.
      const overlay = new THREE.Scene();
      const caustics = new THREE.LineSegments(new THREE.BufferGeometry(), new THREE.LineBasicMaterial({ color: 0xff5050 }));
      const border = new THREE.LineLoop(new THREE.BufferGeometry(), new THREE.LineBasicMaterial({ color: 0x8a93a6 }));
      overlay.add(caustics, border);
      const overlayCamera = new THREE.OrthographicCamera(-1, 1, 1, -1, -1, 1);
      const previousHandler = renderer.debug.onShaderError;
      renderer.debug.onShaderError = (gl, program, vs, fs) => {
        ctx.setNotes('Lensing shader failed to compile:\n' + (gl.getShaderInfoLog(fs) || gl.getProgramInfoLog(program)));
      };

      // The sky (and its spatial index) depends only on the seed and the count: rebuilt only when they change.
      let skyKey = '';
      function buildSky() {
        const key = Math.round(options.sources) + '/' + Math.round(options.seed);
        if (key === skyKey) return;
        skyKey = key;
        const srcs = makeSources(Math.round(options.sources), Math.round(options.seed), 12);
        const n = Math.min(srcs.length, MAXSRC);
        uniforms.uCount.value = n;
        let span = 1;
        for (let i = 0; i < n; i++) {
          const s = srcs[i];
          srcData.set([s.x, s.y, s.size, s.pa, s.q, s.bright, s.type, s.hue, s.pitch, s.arms, s.bulge, s.phase], 12 * i);
          span = Math.max(span, Math.abs(s.x) + 7 * s.size, Math.abs(s.y) + 7 * s.size);
        }
        srcTexture.needsUpdate = true;
        span *= 1.001;
        const cs = 2 * span / GRID, cells = Array.from({ length: GRID * GRID }, () => []);
        const cellOf = (v) => Math.max(0, Math.min(GRID - 1, Math.floor((v + span) / cs)));
        for (let i = 0; i < n; i++) {
          const s = srcs[i], R = 7 * s.size;
          for (let cy = cellOf(s.y - R); cy <= cellOf(s.y + R); cy++) {
            for (let cx = cellOf(s.x - R); cx <= cellOf(s.x + R); cx++) {
              // Keep the cell only if the footprint circle reaches it.
              const nx = Math.max(-span + cx * cs, Math.min(s.x, -span + (cx + 1) * cs)), ny = Math.max(-span + cy * cs, Math.min(s.y, -span + (cy + 1) * cs));
              if ((nx - s.x) ** 2 + (ny - s.y) ** 2 <= R * R) cells[cy * GRID + cx].push(i);
            }
          }
        }
        const list = [];
        cells.forEach((c, k) => { cellData[4 * k] = list.length; cellData[4 * k + 1] = c.length; list.push(...c); });
        cellTexture.needsUpdate = true;
        const rows = Math.max(1, Math.ceil(list.length / IDXW));
        if (!indexTexture || indexTexture.image.height !== rows) {
          if (indexTexture) indexTexture.dispose();
          indexTexture = new THREE.DataTexture(new Float32Array(IDXW * rows * 4), IDXW, rows, THREE.RGBAFormat, THREE.FloatType);
          indexTexture.minFilter = indexTexture.magFilter = THREE.NearestFilter;
          uniforms.uIndex.value = indexTexture;
        }
        const idx = indexTexture.image.data;
        idx.fill(0);
        list.forEach((v, k) => { idx[4 * k] = v; });
        indexTexture.needsUpdate = true;
        uniforms.uGridSpan.value = span;
        skyInfo = `${n} galaxies; at most ${Math.max(...cells.map((c) => c.length))} per pixel (spatial index ${GRID}x${GRID})`;
      }

      let phys = null, lens = null, causticInfo = '', skyInfo = '', lastFrameKey = '', lastChange = 0;
      function update() {
        lastChange = performance.now();
        phys = lensPhysics(options, H0, Om);
        lens = { point: phys.point, b: phys.thetaE, q: options.q, pa: options.pa * Math.PI / 180, gamma: options.gamma, gpa: options.gpa * Math.PI / 180, lx: options._lx || 0, ly: options._ly || 0 };
        // The sky is fixed (it depends only on the seed and count), so galaxies stay put while the lens changes.
        buildSky();
        const stars = makeStars(Math.round(options.seed), 10);
        stars.forEach((st, i) => uniforms.uStar.value[i].set(st.x, st.y, st.bright, 0));
        uniforms.uStars.value = options.stars ? 1 : 0;
        uniforms.uPoint.value = lens.point ? 1 : 0;
        uniforms.uB.value = lens.b; uniforms.uQ.value = lens.q; uniforms.uPA.value = lens.pa;
        uniforms.uGamma.value = lens.gamma; uniforms.uGPA.value = lens.gpa;
        uniforms.uLens.value.set(lens.lx, lens.ly);
        uniforms.uCrit.value = options.crit ? 1 : 0;
        uniforms.uMap.value = Math.round(options.map);
        uniforms.uLensLight.value = options.lensLight;
        uniforms.uExposure.value = options.exposure;
        uniforms.uRing.value = options.ring ? lens.b : 0;
        computeCaustics();
        const notes = [];
        if (options.zs <= options.zl) notes.push('The source must be behind the lens: source redshift raised to lens + 0.01.');
        if (scene.thetaEArcsec) {
          const same = Math.abs(Math.log10(scene.massSolar) - options.mass) < 1e-6 && (scene.model === 'point') === phys.point
            && Math.abs(scene.zLens - options.zl) < 1e-9 && Math.abs(scene.zSource - options.zs) < 1e-9 && Math.abs(scene.apertureKpc - options.aperture) < 1e-9;
          notes.push(same ? `Modelica (GalaxyLensing) gives theta_E = ${scene.thetaEArcsec.toFixed(5)}"; this view: ${phys.thetaE.toFixed(5)}".`
            : `The Modelica scene was computed for 10^${Math.log10(scene.massSolar).toFixed(2)} M_sun, z = ${scene.zLens} -> ${scene.zSource}: theta_E = ${scene.thetaEArcsec.toFixed(4)}".`);
        }
        notes.push('Thin-lens approximation; deflection 4GM/(c^2 b) is the weak-field limit checked by Examples.LightDeflection.');
        ctx.setNotes(notes.join('\n'));
      }
      // Critical curves (det A = 0) by marching squares on a grid, mapped through the lens equation to caustics.
      // They are computed for the lens at the origin and translated with it (the lens equation is
      // translation-invariant), so dragging the lens does not recompute them; the grid size is bucketed in
      // steps of sqrt(2), so zooming recomputes them only now and then.
      let causticKey = '';
      function computeCaustics() {
        caustics.position.set(lens.lx, lens.ly, 0);
        const Rmin = Math.max(3 * lens.b + 1.5 * lens.gamma * 10, 1.5 * options.fov);
        const R = Math.pow(2, Math.ceil(2 * Math.log2(Rmin)) / 2);
        const key = JSON.stringify([lens.point, lens.b, lens.q, lens.pa, lens.gamma, lens.gpa, R]);
        if (key === causticKey) return;
        causticKey = key;
        const L = { ...lens, lx: 0, ly: 0 };
        computeCausticsAt(R, L);
      }
      function computeCausticsAt(R, lens) {
        const N = 220, h = 2 * R / (N - 1);
        const det = new Float64Array(N * N), X = (i) => lens.lx - R + i * h, Y = (j) => lens.ly - R + j * h;
        // One deflection per node of a grid with a one-node border; det A from central differences of
        // neighbouring nodes (a fifth of the deflection evaluations of differencing at each node).
        const M = N + 2, ax = new Float64Array(M * M), ay = new Float64Array(M * M);
        for (let i = 0; i < M; i++) for (let j = 0; j < M; j++) {
          const a = deflection(X(i - 1), Y(j - 1), lens);
          ax[i * M + j] = a[0]; ay[i * M + j] = a[1];
        }
        for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) {
          const e = (i + 2) * M + j + 1, w = i * M + j + 1, n = (i + 1) * M + j + 2, s = (i + 1) * M + j;
          const axx = (ax[e] - ax[w]) / (2 * h), ayx = (ay[e] - ay[w]) / (2 * h), axy = (ax[n] - ax[s]) / (2 * h), ayy = (ay[n] - ay[s]) / (2 * h);
          det[i * N + j] = (1 - axx) * (1 - ayy) - axy * ayx;
        }
        const pts = [];
        const lerp = (p0, p1, f0, f1) => p0 + (p1 - p0) * f0 / (f0 - f1);
        const toSource = (x, y) => { const a = deflection(x, y, lens); return [x - a[0], y - a[1]]; };
        for (let i = 0; i < N - 1; i++) for (let j = 0; j < N - 1; j++) {
          const f00 = det[i * N + j], f10 = det[(i + 1) * N + j], f01 = det[i * N + j + 1], f11 = det[(i + 1) * N + j + 1];
          const cr = [];
          if (f00 * f10 < 0) cr.push([lerp(X(i), X(i + 1), f00, f10), Y(j)]);
          if (f10 * f11 < 0) cr.push([X(i + 1), lerp(Y(j), Y(j + 1), f10, f11)]);
          if (f01 * f11 < 0) cr.push([lerp(X(i), X(i + 1), f01, f11), Y(j + 1)]);
          if (f00 * f01 < 0) cr.push([X(i), lerp(Y(j), Y(j + 1), f00, f01)]);
          // Skip the singular centre, where det changes sign without a critical curve (the cut is added below).
          for (let k = 0; k + 1 < cr.length; k += 2) {
            if (Math.hypot(cr[k][0] - lens.lx, cr[k][1] - lens.ly) < 2.5 * h) continue;
            const s0 = toSource(cr[k][0], cr[k][1]), s1 = toSource(cr[k + 1][0], cr[k + 1][1]);
            pts.push(s0[0], s0[1], 0, s1[0], s1[1], 0);
          }
        }
        // The isothermal galaxy's cut: the image of theta -> 0, beta = lens - alpha(direction), exactly.
        if (!lens.point) {
          const eps = 1e-6 * lens.b, m = 256;
          let prev = null;
          for (let k = 0; k <= m; k++) {
            const a = 2 * Math.PI * k / m, x = lens.lx + eps * Math.cos(a), y = lens.ly + eps * Math.sin(a);
            const cur = toSource(x, y);
            if (prev) pts.push(prev[0], prev[1], 0, cur[0], cur[1], 0);
            prev = cur;
          }
        }
        caustics.geometry.dispose();
        caustics.geometry = new THREE.BufferGeometry();
        caustics.geometry.setAttribute('position', new THREE.Float32BufferAttribute(pts, 3));
        causticInfo = pts.length ? `${pts.length / 6} caustic segments` : 'no caustics in view (weak lensing only)';
      }

      // Drag to move the lens, scroll to zoom.
      const el = renderer.domElement;
      let dragging = null;
      const onDown = (e) => { dragging = { x: e.clientX, y: e.clientY, lx: options._lx || 0, ly: options._ly || 0 }; el.setPointerCapture(e.pointerId); };
      const onMove = (e) => {
        if (!dragging) return;
        const scale = 2 * options.fov / el.clientHeight;
        options._lx = Math.max(-12, Math.min(12, dragging.lx + (e.clientX - dragging.x) * scale));
        options._ly = Math.max(-12, Math.min(12, dragging.ly - (e.clientY - dragging.y) * scale));
        update();
        ctx.refreshReadout();
      };
      const onUp = () => { dragging = null; };
      const onWheel = (e) => {
        e.preventDefault();
        options.fov = Math.max(0.5, Math.min(40, options.fov * Math.exp(0.0015 * e.deltaY)));
        update();
        ctx.refreshReadout();
        clearTimeout(wheelTimer);
        wheelTimer = setTimeout(() => ctx.refreshOptionsUI(), 250);
      };
      let wheelTimer = null;
      el.addEventListener('pointerdown', onDown);
      el.addEventListener('pointermove', onMove);
      el.addEventListener('pointerup', onUp);
      el.addEventListener('wheel', onWheel, { passive: false });

      update();
      return {
        group: new THREE.Group(),
        render(r) {
          // The renderer's logical size (what setViewport expects), not the canvas's CSS box.
          const size = r.getSize(new THREE.Vector2()), w = size.x, h = size.y;
          // The image only changes with the options or the window size: skip unchanged frames.
          // While the lens is being dragged or a slider moves, draw one sample per pixel; 150 ms after the last
          // change, redraw once with 2 x 2 supersampling.
          const settled = performance.now() - lastChange > 150;
          const frameKey = JSON.stringify(options) + w + 'x' + h + r.getPixelRatio() + settled;
          if (frameKey === lastFrameKey) return;
          lastFrameKey = frameKey;
          uniforms.uSuper.value = options.supersample && settled ? 1 : 0;
          const pr = r.getPixelRatio();
          uniforms.uCentre.value.set(0, 0);
          uniforms.uFov.value = options.fov;
          // Main view: the lensed image.
          uniforms.uLensing.value = 1;
          uniforms.uAspect.value = w / h;
          uniforms.uPixel.value = 2 * options.fov / (h * pr);
          r.setScissorTest(false);
          r.setViewport(0, 0, w, h);
          r.render(quadScene, quadCamera);
          if (options.inset) {
            // Inset (bottom right): the unlensed source plane at the same angular scale, with the caustics.
            const s = Math.round(Math.min(w, h) * 0.32), x0 = w - s - 10, y0 = 10;
            r.autoClear = false;
            r.setScissorTest(true);
            r.setScissor(x0, y0, s, s);
            r.setViewport(x0, y0, s, s);
            uniforms.uLensing.value = 0;
            uniforms.uAspect.value = 1;
            uniforms.uPixel.value = 2 * options.fov / (s * pr);
            r.render(quadScene, quadCamera);
            overlayCamera.left = -options.fov; overlayCamera.right = options.fov;
            overlayCamera.top = options.fov; overlayCamera.bottom = -options.fov;
            overlayCamera.updateProjectionMatrix();
            const f = options.fov * 0.999;
            border.geometry.setFromPoints([new THREE.Vector3(-f, -f, 0), new THREE.Vector3(f, -f, 0), new THREE.Vector3(f, f, 0), new THREE.Vector3(-f, f, 0)]);
            r.render(overlay, overlayCamera);
            r.setScissorTest(false);
            r.setViewport(0, 0, w, h);
            r.autoClear = true;
          }
        },
        readout() {
          const lines = [
            `${phys.point ? 'point mass' : 'isothermal galaxy'}: M = 10^${options.mass.toFixed(2)} M_sun${phys.point ? '' : ` within ${options.aperture} kpc, sigma = ${(phys.sigma / 1000).toFixed(0)} km/s`}`,
            `Einstein radius ${phys.thetaE.toFixed(4)}"  (${(phys.Dl * phys.thetaE * ARCSEC / MPC * 1000).toFixed(2)} kpc at the lens)`,
            `mass inside the ring 10^${Math.log10(phys.massInRing).toFixed(2)} M_sun`,
            `D_l ${(phys.Dl / MPC).toFixed(1)}, D_s ${(phys.Ds / MPC).toFixed(1)}, D_ls ${(phys.Dls / MPC).toFixed(1)} Mpc`,
            `Sigma_cr ${phys.SigmaCr.toFixed(2)} kg/m^2;  ${causticInfo}`,
            skyInfo,
          ];
          return lines;
        },
        onOption() { update(); return true; },
        dispose() {
          el.removeEventListener('pointerdown', onDown);
          el.removeEventListener('pointermove', onMove);
          el.removeEventListener('pointerup', onUp);
          el.removeEventListener('wheel', onWheel);
          controls.enabled = true;
          renderer.debug.onShaderError = previousHandler || null;
          material.dispose();
          srcTexture.dispose();
          cellTexture.dispose();
          if (indexTexture) indexTexture.dispose();
          caustics.geometry.dispose();
        },
      };
    },
  });
})();
