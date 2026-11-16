// Ray tracer style: null geodesics integrated per pixel on the GPU, in the spirit of
// https://michaelmoroz.github.io/TracingGeodesics/ : Hamilton's equations for H = 1/2 g^{mu nu} p_mu p_nu,
// with the metric in Kerr-Schild form g = eta + f k k (regular across the horizon, exact inverse
// eta - f K K). Covers the Kerr-Newman family named by the scene's metric record (minkowski,
// schwarzschild, kerr, reissnerNordstrom, kerrNewman); the formula is the one checked in Modelica by
// Metrics.kerrSchild and Tests.GeneralRelativityVisualizationTests.
(function () {
  const FAMILIES = {
    minkowski: () => ({ rs: 0, a: 0, q: 0 }),
    schwarzschild: (p) => ({ rs: p[0], a: 0, q: 0 }),
    kerr: (p) => ({ rs: p[0], a: p[1] || 0, q: 0 }),
    reissnerNordstrom: (p) => ({ rs: p[0], a: 0, q: p[1] || 0 }),
    kerrNewman: (p) => ({ rs: p[0], a: p[1] || 0, q: p[2] || 0 }),
    kerrSchild: (p) => ({ rs: p[0], a: p[1] || 0, q: p[2] || 0 }),
  };

  const VERTEX = /* glsl */`
    varying vec2 vUv;
    void main() { vUv = uv; gl_Position = vec4(position.xy, 0.0, 1.0); }`;

  const FRAGMENT = /* glsl */`
    varying vec2 vUv;
    uniform vec3 uCamPos, uCamRight, uCamUp, uCamFwd;
    uniform float uTanHalfFov, uAspect;
    uniform float uRs, uA, uQ2, uHorizon, uEscape, uStepFrac;
    uniform int uMaxSteps;
    uniform int uDisk, uShift, uSky;
    uniform float uDiskIn, uDiskOut;
    uniform vec2 uJitter;   // sub-pixel offset in uv units, for accumulated anti-aliasing
    uniform int uBodyCount;
    uniform vec4 uBody[4];
    uniform vec3 uBodyColor[4];
    const float PI = 3.14159265358979;

    float blR(vec3 x) {
      float a2 = uA * uA;
      float w = dot(x, x) - a2;
      return sqrt(0.5 * (w + sqrt(w * w + 4.0 * a2 * x.z * x.z)));
    }
    // Kerr-Schild data: covariant null covector k and f = r^2 (rs r - rQ^2) / (r^4 + a^2 z^2).
    void ks(vec3 x, out float f, out vec4 k) {
      float r = max(blR(x), 1e-5);
      float r2 = r * r;
      float d = r2 + uA * uA;
      k = vec4(1.0, (r * x.x + uA * x.y) / d, (r * x.y - uA * x.x) / d, x.z / r);
      f = r2 * (uRs * r - uQ2) / (r2 * r2 + uA * uA * x.z * x.z);
    }
    // Kerr-Schild data with analytic spatial gradients: f, grad f and grad k_i (rows: d/dx, d/dy, d/dz).
    // One evaluation replaces the six of a central difference; checked against finite differences.
    void ksGrad(vec3 x, out float f, out vec4 k, out vec3 gf, out vec3 gk1, out vec3 gk2, out vec3 gk3) {
      float a = uA, a2 = a * a;
      float w = dot(x, x) - a2;
      float D = sqrt(w * w + 4.0 * a2 * x.z * x.z);
      float r2 = max(0.5 * (w + D), 1e-10);
      float r = sqrt(r2);
      vec3 gw = 2.0 * x;
      vec3 gD = (w * gw + vec3(0.0, 0.0, 4.0 * a2 * x.z)) / max(D, 1e-12);
      vec3 gr2 = 0.5 * (gw + gD);
      vec3 gr = gr2 / (2.0 * r);
      float d = r2 + a2;
      float n1 = r * x.x + a * x.y, n2 = r * x.y - a * x.x;
      k = vec4(1.0, n1 / d, n2 / d, x.z / r);
      gk1 = ((gr * x.x + vec3(r, a, 0.0)) * d - n1 * gr2) / (d * d);
      gk2 = ((gr * x.y + vec3(-a, r, 0.0)) * d - n2 * gr2) / (d * d);
      gk3 = (vec3(0.0, 0.0, r) - x.z * gr) / r2;
      float N = r2 * (uRs * r - uQ2), Q = r2 * r2 + a2 * x.z * x.z;
      vec3 gN = (3.0 * uRs * r2 - 2.0 * uQ2 * r) * gr;
      vec3 gQ = 4.0 * r2 * r * gr + vec3(0.0, 0.0, 2.0 * a2 * x.z);
      f = N / Q;
      gf = (gN * Q - N * gQ) / (Q * Q);
    }
    // Hamilton's equations for H = 1/2 g^{mu nu} p_mu p_nu, g^{-1} = eta - f K K, K^mu = (-k_0, k_i):
    // vel = dx/d lambda = dH/dp, frc = dp/d lambda = -dH/dx.
    void derivs(vec3 x, vec4 p, out vec3 vel, out vec3 frc) {
      float f; vec4 k; vec3 gf, gk1, gk2, gk3;
      ksGrad(x, f, k, gf, gk1, gk2, gk3);
      float kp = -k.x * p.x + dot(k.yzw, p.yzw);
      vel = p.yzw - f * kp * k.yzw;
      frc = 0.5 * gf * kp * kp + f * kp * (gk1 * p.y + gk2 * p.z + gk3 * p.w);
    }
    vec3 velocity(vec3 x, vec4 p) {
      vec3 vel, frc;
      derivs(x, p, vel, frc);
      return vel;
    }

    float hash13(vec3 p) {
      p = fract(p * 0.1031);
      p += dot(p, p.zyx + 31.32);
      return fract((p.x + p.y) * p.z);
    }
    vec3 starLayer(vec3 dir, float cells, float density) {
      vec3 q = dir * cells;
      vec3 cell = floor(q);
      float h = hash13(cell);
      if (h < 1.0 - density) return vec3(0.0);
      vec3 centre = cell + 0.5 + 0.35 * (vec3(hash13(cell + 1.7), hash13(cell + 3.1), hash13(cell + 5.3)) - 0.5);
      float d = length(q - centre);
      float b = (h - (1.0 - density)) / density;
      vec3 tint = mix(vec3(1.0, 0.8, 0.6), vec3(0.7, 0.8, 1.0), hash13(cell + 9.1));
      return tint * b * 2.5 * smoothstep(0.45, 0.0, d);
    }
    vec3 sky(vec3 dir) {
      float th = acos(clamp(dir.z, -1.0, 1.0));
      float ph = atan(dir.y, dir.x);
      vec3 col = vec3(0.0006, 0.0008, 0.0016);
      float band = dot(dir, normalize(vec3(0.25, 0.55, 0.8)));
      col += vec3(0.010, 0.009, 0.014) * exp(-10.0 * band * band);
      col += starLayer(dir, 90.0, 0.03) + starLayer(dir, 220.0, 0.02) * 0.6;
      if (uSky == 1) {
        float spacing = PI / 12.0;
        float lat = abs(fract(th / spacing + 0.5) - 0.5) * spacing;
        float lon = abs(fract(ph / spacing + 0.5) - 0.5) * spacing * sin(th);
        float line = min(lat, lon);
        col += vec3(0.03, 0.07, 0.12) * smoothstep(0.004, 0.0, line);
      } else if (uSky == 2) {
        float c = mod(floor(th / (PI / 12.0)) + floor((ph + PI) / (PI / 12.0)), 2.0);
        col = mix(vec3(0.01, 0.012, 0.025), vec3(0.08, 0.09, 0.13), c) + col * 0.5;
      }
      return col;
    }
    // Approximate blackbody colour (Tanner Helland's fit), temperature in kelvin.
    vec3 blackbody(float T) {
      T = clamp(T, 1000.0, 40000.0) / 100.0;
      vec3 c;
      c.r = T <= 66.0 ? 1.0 : clamp(1.292936 * pow(T - 60.0, -0.1332047), 0.0, 1.0);
      c.g = T <= 66.0 ? clamp(0.3900816 * log(T) - 0.6318414, 0.0, 1.0) : clamp(1.1298909 * pow(T - 60.0, -0.0755148), 0.0, 1.0);
      c.b = T >= 66.0 ? 1.0 : (T <= 19.0 ? 0.0 : clamp(0.5432068 * log(T - 10.0) - 1.1962541, 0.0, 1.0));
      return c;
    }
    vec3 diskColour(vec3 h, vec4 p, float puObs) {
      float r = blR(h);
      float x = uDiskIn / r;
      // Thin-disk temperature profile r^-3/4 (1 - sqrt(r_in/r))^1/4, peak normalized to 1.
      float temp = pow(x, 0.75) * pow(max(1.0 - sqrt(x), 0.0), 0.25) / 0.488;
      float g = 1.0;
      if (uShift == 1) {
        float M = 0.5 * uRs;
        float Mr = max(M * r - uQ2, 0.0);
        float omega = sqrt(Mr) / (r * r + uA * sqrt(Mr));
        vec4 v = vec4(1.0, -omega * h.y, omega * h.x, 0.0);
        float f; vec4 k; ks(h, f, k);
        float kv = dot(k, v);
        float gvv = -1.0 + dot(v.yzw, v.yzw) + f * kv * kv;
        if (gvv < 0.0) {
          vec4 u = v / sqrt(-gvv);
          g = puObs / dot(p, u);
        }
      }
      vec3 c = blackbody(4200.0 * g * temp);
      return c * pow(g, 4.0) * pow(temp, 4.0) * 0.9;
    }

    void main() {
      vec2 uv = vUv * 2.0 - 1.0 + uJitter;
      vec3 d = normalize(uCamFwd + uv.x * uTanHalfFov * uAspect * uCamRight + uv.y * uTanHalfFov * uCamUp);
      // Null momentum at the camera: contravariant p = (P0, d), past-directed so the ray runs back along
      // the light that reached the camera; then lower it with g = eta + f k k.
      float fc; vec4 kc; ks(uCamPos, fc, kc);
      float s = dot(kc.yzw, d);
      float A = fc - 1.0, B = 2.0 * fc * s, C = 1.0 + fc * s * s;
      float P0 = (-B + sqrt(max(B * B - 4.0 * A * C, 0.0))) / (2.0 * A);
      vec4 p = vec4(-P0, d) + fc * (P0 + s) * kc;
      float puObs = p.x / sqrt(max(1.0 - fc, 1e-6));   // static observer at the camera

      vec3 x = uCamPos;
      vec3 colour = vec3(0.0);
      float transmit = 1.0;
      bool done = false;
      vec3 v = velocity(x, p);
      for (int i = 0; i < 4096; i++) {
        if (i >= uMaxSteps) break;
        float r = blR(x);
        if (uHorizon > 0.0 && r < uHorizon * 1.01) { done = true; break; }
        if (r > uEscape && dot(x, v) > 0.0) { colour += transmit * sky(normalize(v)); done = true; break; }
        float ds = uStepFrac * max(r - uHorizon, 0.02 * r + 1e-4);
        float dl = ds / max(length(v), 1e-8);
        vec3 k1x, k1p, k2x, k2p, k3x, k3p, k4x, k4p;
        derivs(x, p, k1x, k1p);
        vec3 x2 = x + 0.5 * dl * k1x;    vec4 p2 = vec4(p.x, p.yzw + 0.5 * dl * k1p);
        derivs(x2, p2, k2x, k2p);
        vec3 x3 = x + 0.5 * dl * k2x;    vec4 p3 = vec4(p.x, p.yzw + 0.5 * dl * k2p);
        derivs(x3, p3, k3x, k3p);
        vec3 x4 = x + dl * k3x;          vec4 p4 = vec4(p.x, p.yzw + dl * k3p);
        derivs(x4, p4, k4x, k4p);
        vec3 xn = x + dl / 6.0 * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
        vec4 pn = vec4(p.x, p.yzw + dl / 6.0 * (k1p + 2.0 * k2p + 2.0 * k3p + k4p));
        // Bodies: straight segment x -> xn against each sphere.
        vec3 seg = xn - x;
        float segLen = length(seg);
        vec3 sd = seg / max(segLen, 1e-12);
        float hitT = 1e30; int hitBody = -1;
        for (int b = 0; b < 4; b++) {
          if (b >= uBodyCount) break;
          vec3 oc = x - uBody[b].xyz;
          float bq = dot(oc, sd);
          float cq = dot(oc, oc) - uBody[b].w * uBody[b].w;
          float disc = bq * bq - cq;
          if (disc > 0.0) {
            float t0 = -bq - sqrt(disc);
            if (t0 >= 0.0 && t0 <= segLen && t0 < hitT) { hitT = t0; hitBody = b; }
          }
        }
        // Thin disk in the equatorial plane z = 0.
        float diskT = 1e30;
        if (uDisk == 1 && x.z * xn.z < 0.0) {
          float t = x.z / (x.z - xn.z);
          vec3 hp = mix(x, xn, t);
          float rh = blR(hp);
          if (rh >= uDiskIn && rh <= uDiskOut) diskT = t * segLen;
        }
        if (diskT < hitT) {
          float t = diskT / max(segLen, 1e-12);
          colour += transmit * diskColour(mix(x, xn, t), mix(p, pn, t), puObs);
          transmit *= 0.12;
        }
        if (hitBody >= 0) {
          vec3 hp = x + sd * hitT;
          vec3 n = normalize(hp - uBody[hitBody].xyz);
          float lit = 0.35 + 0.65 * max(dot(n, -sd), 0.0);
          colour += transmit * uBodyColor[hitBody] * lit * 1.4;
          done = true;
          break;
        }
        x = xn; p = pn;
        v = velocity(x, p);
        if (transmit < 0.02) { done = true; break; }
      }
      // Out of steps: a ray still near the hole is winding onto the horizon (traced backward in time it
      // approaches the past horizon only asymptotically in ingoing Kerr-Schild coordinates), so it is
      // captured; a ray still far out is given the sky in its current direction.
      if (!done && blR(x) > 4.0 * max(uHorizon, 1e-6)) colour += transmit * sky(normalize(v));
      // Linear radiance with alpha 1: samples are summed in a float buffer and averaged by alpha at display.
      gl_FragColor = vec4(colour, 1.0);
    }`;

  window.GRStyles = window.GRStyles || [];
  window.GRStyles.push({
    id: 'raytrace',
    label: 'Ray tracer',
    order: 5,
    description: 'Null geodesics traced per pixel on the GPU: the shadow and photon ring, the lensed sky, a thin accretion disk with Doppler and gravitational shifts, and the bodies as a camera sees them. Drag to orbit; or ride along with a body.',
    options: [
      { id: 'resolution', label: 'max resolution', type: 'range', min: 0.25, max: 1, step: 0.05, value: 1, help: 'Fraction of the display\'s device pixels traced once the view is still' },
      { id: 'samples', label: 'refine samples', type: 'range', min: 1, max: 64, step: 1, value: 24, help: 'Jittered full-resolution samples averaged while the view is still (anti-aliasing)' },
      { id: 'steps', label: 'max steps', type: 'range', min: 100, max: 2000, step: 50, value: 600 },
      { id: 'accuracy', label: 'step size', type: 'range', min: 0.01, max: 0.12, step: 0.005, value: 0.04, help: 'Spatial step as a fraction of the distance to the horizon' },
      { id: 'disk', label: 'accretion disk', type: 'checkbox', value: 1 },
      { id: 'diskOuter', label: 'disk outer (rs)', type: 'range', min: 3, max: 30, step: 0.5, value: 8 },
      { id: 'shift', label: 'Doppler + redshift', type: 'checkbox', value: 1 },
      { id: 'sky', label: 'sky', type: 'select', choices: [{ value: 0, label: 'stars' }, { value: 1, label: 'stars + grid' }, { value: 2, label: 'checker' }], value: 1 },
      { id: 'bodySize', label: 'body size (rs)', type: 'range', min: 0.05, max: 1.5, step: 0.05, value: 0.35 },
      { id: 'camera', label: 'camera', type: 'select', value: 0,
        choices: (data) => [{ value: 0, label: 'orbit (drag)' }].concat(data.trackSpecs.map((t, i) => ({ value: i + 1, label: 'ride with ' + (t.name || 'body ' + (i + 1)) }))) },
      { id: 'paths', label: 'overlay paths', type: 'checkbox', value: 0, help: 'Draw the trajectories on top, without lensing' },
      { id: 'exposure', label: 'exposure', type: 'range', min: 0.2, max: 4, step: 0.1, value: 1.0 },
    ],
    available(data) {
      if (!data.metric) return 'needs a metric record (Visualization.writeMetric, or SceneRecorder metricName)';
      if (!FAMILIES[data.metric.name]) return `the ray tracer does not know the metric "${data.metric.name}"`;
      return null;
    },
    create(ctx) {
      const { THREE, data, options, renderer, camera, controls } = ctx;
      const params = FAMILIES[data.metric.name](data.metric.p || []);
      const M = 0.5 * params.rs;
      const disc = M * M - params.a * params.a - params.q * params.q;
      const horizon = params.rs > 0 && disc >= 0 ? M + Math.sqrt(disc) : 0;
      // Prograde ISCO of Kerr (Bardeen, Press and Teukolsky); a charge is ignored here.
      let isco = 0;
      if (M > 0) {
        const chi = Math.min(Math.abs(params.a) / M, 0.999);
        const z1 = 1 + Math.cbrt(1 - chi * chi) * (Math.cbrt(1 + chi) + Math.cbrt(1 - chi));
        const z2 = Math.sqrt(3 * chi * chi + z1 * z1);
        isco = M * (3 + z2 - Math.sqrt((3 - z1) * (3 + z1 + 2 * z2)));
      }
      const unit = params.rs > 0 ? params.rs : ctx.referenceRadius;
      const notes = [`${data.metric.name}: rs = ${ctx.fmt(params.rs)}, a = ${ctx.fmt(params.a)}, rQ = ${ctx.fmt(params.q)}; horizon r+ = ${ctx.fmt(horizon)}, disk from the ISCO r = ${ctx.fmt(isco)}.`,
        'Snapshot of a stationary metric: bodies are drawn where they are now (no light-travel delay); the camera is a static observer.'];
      if (params.a !== 0) notes.push('Boyer-Lindquist positions are mapped to Kerr-Schild x + iy = (r + ia) sin(theta) e^(i phi), ignoring the small azimuth shift.');
      ctx.setNotes(notes.join('\n'));

      const uniforms = {
        uCamPos: { value: new THREE.Vector3() }, uCamRight: { value: new THREE.Vector3() },
        uCamUp: { value: new THREE.Vector3() }, uCamFwd: { value: new THREE.Vector3() },
        uTanHalfFov: { value: 0.5 }, uAspect: { value: 1 },
        uRs: { value: params.rs }, uA: { value: params.a }, uQ2: { value: params.q * params.q },
        uHorizon: { value: horizon }, uEscape: { value: 100 }, uStepFrac: { value: options.accuracy },
        uMaxSteps: { value: options.steps }, uDisk: { value: 1 }, uShift: { value: 1 }, uSky: { value: 1 },
        uDiskIn: { value: isco || unit }, uDiskOut: { value: 8 * unit }, uJitter: { value: new THREE.Vector2() },
        uBodyCount: { value: 0 },
        uBody: { value: [new THREE.Vector4(), new THREE.Vector4(), new THREE.Vector4(), new THREE.Vector4()] },
        uBodyColor: { value: [new THREE.Color(), new THREE.Color(), new THREE.Color(), new THREE.Color()] },
      };
      // Two materials share the shader and uniforms: one overwrites (previews), one adds into the
      // accumulation buffer (refinement samples), counting samples in alpha.
      const material = new THREE.ShaderMaterial({ vertexShader: VERTEX, fragmentShader: FRAGMENT, uniforms, depthTest: false, depthWrite: false, blending: THREE.NoBlending });
      const addMaterial = new THREE.ShaderMaterial({ vertexShader: VERTEX, fragmentShader: FRAGMENT, uniforms, depthTest: false, depthWrite: false,
        transparent: true, blending: THREE.CustomBlending, blendEquation: THREE.AddEquation, blendSrc: THREE.OneFactor, blendDst: THREE.OneFactor,
        blendEquationAlpha: THREE.AddEquation, blendSrcAlpha: THREE.OneFactor, blendDstAlpha: THREE.OneFactor });
      const quad = new THREE.Mesh(new THREE.PlaneGeometry(2, 2), material);
      const quadScene = new THREE.Scene();
      quadScene.add(quad);
      const quadCamera = new THREE.OrthographicCamera(-1, 1, 1, -1, 0, 1);
      // Display: average the accumulated samples (rgb / alpha), fall back to the preview where a band is
      // not traced yet, then tone-map with the exposure (so exposure never forces a re-trace).
      const displayMaterial = new THREE.ShaderMaterial({
        vertexShader: VERTEX,
        fragmentShader: `varying vec2 vUv;
          uniform sampler2D tAccum, tPreview; uniform float uUseAccum, uExposure;
          void main() {
            vec4 a = texture2D(tAccum, vUv);
            vec4 pv = texture2D(tPreview, vUv);
            vec3 c = (uUseAccum > 0.5 && a.a > 0.5) ? a.rgb / a.a : pv.rgb / max(pv.a, 1e-6);
            c = vec3(1.0) - exp(-c * uExposure);
            gl_FragColor = vec4(pow(c, vec3(1.0 / 2.2)), 1.0);
          }`,
        uniforms: { tAccum: { value: null }, tPreview: { value: null }, uUseAccum: { value: 0 }, uExposure: { value: 1 } },
        depthTest: false, depthWrite: false,
      });
      const blitScene = new THREE.Scene();
      blitScene.add(new THREE.Mesh(new THREE.PlaneGeometry(2, 2), displayMaterial));
      // Half-float targets hold linear radiance and sums of samples; without float rendering only previews are used.
      const floatOk = renderer.capabilities.isWebGL2 && renderer.extensions.has('EXT_color_buffer_float');
      const texType = floatOk ? THREE.HalfFloatType : THREE.UnsignedByteType;
      const makeTarget = (w, h) => new THREE.WebGLRenderTarget(w, h, { type: texType, minFilter: THREE.LinearFilter, magFilter: THREE.LinearFilter, depthBuffer: false });
      let preview = null, accum = null;
      // Refinement state.
      let lastKey = '', lastChange = 0, level = 0, sample = 0, band = 0, bands = 1;
      let pixelBudget = 350000, lastFrame = performance.now(), status = '';
      const halton = (i, b) => { let f = 1, r = 0; while (i > 0) { f /= b; r += f * (i % b); i = Math.floor(i / b); } return r; };

      // Camera damping keeps easing the camera for seconds after the mouse is released, which would keep the
      // tracer in preview mode; the ray tracer stops the camera when the mouse stops.
      const previousDamping = controls.enableDamping;
      controls.enableDamping = false;

      // Shader compile errors are shown in the side panel rather than only in the console.
      const previousHandler = renderer.debug.onShaderError;
      renderer.debug.onShaderError = (gl, program, vs, fs) => {
        const log = gl.getShaderInfoLog(fs) || gl.getShaderInfoLog(vs) || gl.getProgramInfoLog(program);
        ctx.setNotes('Ray-tracer shader failed to compile:\n' + log);
      };

      // World (three.js, y up) <-> physical Kerr-Schild Cartesian (z = spin axis).
      const toPhys = (w) => new THREE.Vector3(w.x, -w.z, w.y);
      const samplePhys = (x) => {
        if (data.header.chart === 'spherical') {
          const rho = Math.sqrt(x[1] * x[1] + params.a * params.a) * Math.sin(x[2]);
          return new THREE.Vector3(rho * Math.cos(x[3]), rho * Math.sin(x[3]), x[1] * Math.cos(x[2]));
        }
        return new THREE.Vector3(x[1], x[2], x[3]);
      };
      const physToWorld = (p) => new THREE.Vector3(p.x, p.z, -p.y);

      function syncUniforms() {
        const bodies = ctx.bodiesNow();
        if (options.camera > 0 && bodies[options.camera - 1]) {
          // Ride along: just outside the chosen body, looking at the hole.
          const P = samplePhys(bodies[options.camera - 1].sample.x);
          const out = P.clone().normalize().multiplyScalar(3 * options.bodySize * unit);
          camera.position.copy(physToWorld(P.clone().add(out)));
          camera.up.set(0, 1, 0);
          camera.lookAt(0, 0, 0);
          controls.enabled = false;
        } else {
          controls.enabled = true;
        }
        camera.updateMatrixWorld();
        const e = camera.matrixWorld.elements;
        uniforms.uCamPos.value.copy(toPhys(camera.position));
        uniforms.uCamRight.value.copy(toPhys(new THREE.Vector3(e[0], e[1], e[2]).normalize()));
        uniforms.uCamUp.value.copy(toPhys(new THREE.Vector3(e[4], e[5], e[6]).normalize()));
        uniforms.uCamFwd.value.copy(toPhys(new THREE.Vector3(-e[8], -e[9], -e[10]).normalize()));
        uniforms.uTanHalfFov.value = Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2);
        uniforms.uAspect.value = camera.aspect;
        uniforms.uEscape.value = Math.max(2 * uniforms.uCamPos.value.length(), 40 * unit, 2 * ctx.extent);
        uniforms.uStepFrac.value = options.accuracy;
        uniforms.uMaxSteps.value = Math.round(options.steps);
        uniforms.uDisk.value = options.disk ? 1 : 0;
        uniforms.uDiskOut.value = options.diskOuter * unit;
        uniforms.uShift.value = options.shift ? 1 : 0;
        uniforms.uSky.value = Math.round(options.sky);
        // Bodies hidden in the legend are not traced.
        const shown = bodies.filter(b => !b.hidden).slice(0, 4);
        const count = shown.length;
        uniforms.uBodyCount.value = count;
        for (let i = 0; i < count; i++) {
          const P = samplePhys(shown[i].sample.x);
          uniforms.uBody.value[i].set(P.x, P.y, P.z, options.bodySize * unit);
          uniforms.uBodyColor.value[i].copy(shown[i].color);
        }
        // Everything that changes the image, to skip tracing when nothing moved.
        const { exposure, paths, ...traced } = options;
        // Rounded so that motion far below a pixel (the tail of an easing camera) does not restart refinement.
        return [e.map(v => v.toPrecision(5)).join(','), uniforms.uBody.value.slice(0, count).map(b => b.toArray().join(',')).join(';'),
          JSON.stringify(traced), renderer.domElement.width, renderer.domElement.height].join('|');
      }

      // Size of a pass: previews are fractions of CSS pixels, refinement uses device pixels.
      function passSize(levelIndex) {
        const el = renderer.domElement;
        if (levelIndex === 0) return [Math.max(16, Math.round(el.clientWidth * 0.3)), Math.max(16, Math.round(el.clientHeight * 0.3))];
        if (levelIndex === 1) return [Math.max(16, Math.round(el.clientWidth * 0.6)), Math.max(16, Math.round(el.clientHeight * 0.6))];
        return [Math.max(16, Math.round(el.width * options.resolution)), Math.max(16, Math.round(el.height * options.resolution))];
      }
      function ensure(t, w, h) {
        if (t && t.width === w && t.height === h) return t;
        if (t) t.dispose();
        return makeTarget(w, h);
      }
      function trace(r, targetRT, mat, stepFrac, maxSteps, jx, jy, y0, h) {
        uniforms.uStepFrac.value = stepFrac;
        uniforms.uMaxSteps.value = maxSteps;
        uniforms.uJitter.value.set(jx, jy);
        quad.material = mat;
        if (h !== undefined) { targetRT.scissor.set(0, y0, targetRT.width, h); targetRT.scissorTest = true; }
        else targetRT.scissorTest = false;
        r.setRenderTarget(targetRT);
        r.render(quadScene, quadCamera);
        r.setRenderTarget(null);
        targetRT.scissorTest = false;
      }
      // One frame of progressive refinement: a coarse preview while anything changes, a sharper preview once
      // still, then full-resolution samples with finer steps and sub-pixel jitter, traced in horizontal bands
      // that fit a per-frame pixel budget and summed into the accumulation buffer.
      function refine(r, key) {
        const now = performance.now();
        const frameTime = now - lastFrame;
        lastFrame = now;
        if (key !== lastKey) { lastKey = key; lastChange = now; level = 0; }
        const still = now - lastChange > 150;
        if (level === 0 || !still) {
          const [w, h] = passSize(0);
          preview = ensure(preview, w, h);
          trace(r, preview, material, options.accuracy * 1.5, Math.round(options.steps), 0, 0);
          if (still) level = 1;
          status = `preview ${w} x ${h}`;
          return false;
        }
        if (level === 1) {
          const [w, h] = passSize(1);
          preview = ensure(preview, w, h);
          trace(r, preview, material, options.accuracy, Math.round(options.steps), 0, 0);
          level = floatOk ? 2 : 3;
          sample = 0; band = 0;
          status = `preview ${w} x ${h}`;
          return false;
        }
        if (level === 3) return sample > 0;   // converged (or no float targets)
        const [w, h] = passSize(2);
        if (!accum || accum.width !== w || accum.height !== h || (sample === 0 && band === 0)) {
          accum = ensure(accum, w, h);
          const old = r.getClearAlpha(), oldColour = r.getClearColor(new THREE.Color());
          r.setRenderTarget(accum);
          r.setClearColor(0x000000, 0);
          r.clear(true, false, false);
          r.setRenderTarget(null);
          r.setClearColor(oldColour, old);
          // Adapt the per-frame budget to how long the previous refinement frame took.
          if (frameTime > 70) pixelBudget = Math.max(60000, pixelBudget * 0.7);
          bands = Math.max(1, Math.ceil(w * h / pixelBudget));
        } else if (band === 0) {
          if (frameTime > 70) pixelBudget = Math.max(60000, pixelBudget * 0.7);
          else if (frameTime < 30) pixelBudget = Math.min(4e6, pixelBudget * 1.25);
          bands = Math.max(1, Math.ceil(w * h / pixelBudget));
        }
        const bandH = Math.ceil(h / bands), y0 = band * bandH;
        // Halton (2, 3) jitter within the pixel; the first sample is centred.
        const jx = sample === 0 ? 0 : (halton(sample, 2) - 0.5) * 2 / w;
        const jy = sample === 0 ? 0 : (halton(sample, 3) - 0.5) * 2 / h;
        trace(r, accum, addMaterial, options.accuracy * 0.6, Math.round(options.steps * 1.5), jx, jy, y0, Math.min(bandH, h - y0));
        band++;
        if (band >= bands) { band = 0; sample++; }
        if (sample >= options.samples) level = 3;
        status = level === 3 ? `refined: ${sample} samples at ${w} x ${h}` : `refining: sample ${sample + 1} of ${options.samples} at ${w} x ${h}`;
        return true;
      }

      return {
        group: new THREE.Group(),
        hidesBodies: true,
        place(x) { return physToWorld(samplePhys(x)); },
        camera() {
          const R = Math.min(Math.max(14 * unit, 0.9 * ctx.extent), 60 * unit);
          // Slightly above the disk plane, so the disk's far side is lensed over the shadow.
          const incl = THREE.MathUtils.degToRad(80);
          return { position: new THREE.Vector3(R * Math.sin(incl), R * Math.cos(incl), R * 0.15), target: new THREE.Vector3(0, 0, 0) };
        },
        render(r, cam, scene3) {
          const key = syncUniforms();
          const useAccum = refine(r, key);
          displayMaterial.uniforms.tAccum.value = accum ? accum.texture : preview.texture;
          displayMaterial.uniforms.tPreview.value = preview.texture;
          displayMaterial.uniforms.uUseAccum.value = useAccum && accum ? 1 : 0;
          displayMaterial.uniforms.uExposure.value = options.exposure;
          r.render(blitScene, quadCamera);
          if (options.paths) {
            const background = scene3.background;
            scene3.background = null;
            const autoClear = r.autoClear;
            r.autoClear = false;
            r.clearDepth();
            r.render(scene3, cam);
            r.autoClear = autoClear;
            scene3.background = background;
          }
        },
        onOption(id) {
          if (id === 'camera') { controls.enabled = !(options.camera > 0); return true; }
          return true;   // every other option is a uniform, picked up on the next frame
        },
        readout() {
          return [`horizon r+ ${ctx.fmt(horizon)}, ISCO ${ctx.fmt(isco)}`, status];
        },
        dispose() {
          controls.enabled = true;
          controls.enableDamping = previousDamping;
          renderer.debug.onShaderError = previousHandler || null;
          if (preview) preview.dispose();
          if (accum) accum.dispose();
          material.dispose();
          addMaterial.dispose();
          displayMaterial.dispose();
        },
      };
    },
  });
})();
