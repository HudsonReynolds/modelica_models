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
    uniform float uDiskIn, uDiskOut, uExposure;
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
    // H = 1/2 g^{mu nu} p_mu p_nu with g^{-1} = eta - f K K, K^mu = (-k_0, k_i).
    float hamiltonian(vec3 x, vec4 p) {
      float f; vec4 k; ks(x, f, k);
      float kp = -k.x * p.x + dot(k.yzw, p.yzw);
      return 0.5 * (-p.x * p.x + dot(p.yzw, p.yzw) - f * kp * kp);
    }
    vec3 velocity(vec3 x, vec4 p) {
      float f; vec4 k; ks(x, f, k);
      float kp = -k.x * p.x + dot(k.yzw, p.yzw);
      return p.yzw - f * kp * k.yzw;
    }
    vec3 force(vec3 x, vec4 p) {
      float h = 2e-3 * max(blR(x), 0.25 * max(uRs, 1e-3));
      vec3 g;
      g.x = hamiltonian(x + vec3(h, 0.0, 0.0), p) - hamiltonian(x - vec3(h, 0.0, 0.0), p);
      g.y = hamiltonian(x + vec3(0.0, h, 0.0), p) - hamiltonian(x - vec3(0.0, h, 0.0), p);
      g.z = hamiltonian(x + vec3(0.0, 0.0, h), p) - hamiltonian(x - vec3(0.0, 0.0, h), p);
      return -g / (2.0 * h);
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
      vec2 uv = vUv * 2.0 - 1.0;
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
        vec3 k1x = v;                    vec3 k1p = force(x, p);
        vec3 x2 = x + 0.5 * dl * k1x;    vec4 p2 = vec4(p.x, p.yzw + 0.5 * dl * k1p);
        vec3 k2x = velocity(x2, p2);     vec3 k2p = force(x2, p2);
        vec3 x3 = x + 0.5 * dl * k2x;    vec4 p3 = vec4(p.x, p.yzw + 0.5 * dl * k2p);
        vec3 k3x = velocity(x3, p3);     vec3 k3p = force(x3, p3);
        vec3 x4 = x + dl * k3x;          vec4 p4 = vec4(p.x, p.yzw + dl * k3p);
        vec3 k4x = velocity(x4, p4);     vec3 k4p = force(x4, p4);
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
      colour = vec3(1.0) - exp(-colour * uExposure);
      gl_FragColor = vec4(pow(colour, vec3(1.0 / 2.2)), 1.0);
    }`;

  window.GRStyles = window.GRStyles || [];
  window.GRStyles.push({
    id: 'raytrace',
    label: 'Ray tracer',
    order: 5,
    description: 'Null geodesics traced per pixel on the GPU: the shadow and photon ring, the lensed sky, a thin accretion disk with Doppler and gravitational shifts, and the bodies as a camera sees them. Drag to orbit; or ride along with a body.',
    options: [
      { id: 'resolution', label: 'resolution', type: 'range', min: 0.2, max: 1, step: 0.05, value: 0.6, help: 'Fraction of the window resolution traced when the view is still' },
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
        uDiskIn: { value: isco || unit }, uDiskOut: { value: 8 * unit }, uExposure: { value: 1 },
        uBodyCount: { value: 0 },
        uBody: { value: [new THREE.Vector4(), new THREE.Vector4(), new THREE.Vector4(), new THREE.Vector4()] },
        uBodyColor: { value: [new THREE.Color(), new THREE.Color(), new THREE.Color(), new THREE.Color()] },
      };
      const material = new THREE.ShaderMaterial({ vertexShader: VERTEX, fragmentShader: FRAGMENT, uniforms, depthTest: false, depthWrite: false });
      const quadScene = new THREE.Scene();
      const quadCamera = new THREE.OrthographicCamera(-1, 1, 1, -1, 0, 1);
      quadScene.add(new THREE.Mesh(new THREE.PlaneGeometry(2, 2), material));
      const blitMaterial = new THREE.ShaderMaterial({
        vertexShader: VERTEX,
        fragmentShader: 'varying vec2 vUv; uniform sampler2D tImage; void main() { gl_FragColor = texture2D(tImage, vUv); }',
        uniforms: { tImage: { value: null } }, depthTest: false, depthWrite: false,
      });
      const blitScene = new THREE.Scene();
      blitScene.add(new THREE.Mesh(new THREE.PlaneGeometry(2, 2), blitMaterial));
      let target = null;
      let lastKey = '', lastChange = 0, renderedScale = 0;

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
        uniforms.uExposure.value = options.exposure;
        const count = Math.min(bodies.length, 4);
        uniforms.uBodyCount.value = count;
        for (let i = 0; i < count; i++) {
          const P = samplePhys(bodies[i].sample.x);
          uniforms.uBody.value[i].set(P.x, P.y, P.z, options.bodySize * unit);
          uniforms.uBodyColor.value[i].copy(bodies[i].color);
        }
        // Everything that changes the image, to skip tracing when nothing moved.
        return [e.join(','), uniforms.uBody.value.slice(0, count).map(b => b.toArray().join(',')).join(';'),
          JSON.stringify(options), renderer.domElement.width, renderer.domElement.height].join('|');
      }

      function ensureTarget(scale) {
        const w = Math.max(16, Math.round(renderer.domElement.clientWidth * scale));
        const h = Math.max(16, Math.round(renderer.domElement.clientHeight * scale));
        if (!target || target.width !== w || target.height !== h) {
          if (target) target.dispose();
          target = new THREE.WebGLRenderTarget(w, h, { minFilter: THREE.LinearFilter, magFilter: THREE.LinearFilter, depthBuffer: false });
        }
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
          const now = performance.now();
          if (key !== lastKey) { lastKey = key; lastChange = now; }
          const moving = now - lastChange < 250;
          const scale = moving ? Math.min(options.resolution, 0.35) : options.resolution;
          if (key !== material.userData.key || scale !== renderedScale) {
            ensureTarget(scale);
            r.setRenderTarget(target);
            r.render(quadScene, quadCamera);
            r.setRenderTarget(null);
            material.userData.key = key;
            renderedScale = scale;
          }
          blitMaterial.uniforms.tImage.value = target.texture;
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
          return [`horizon r+ ${ctx.fmt(horizon)}, ISCO ${ctx.fmt(isco)}, traced at ${Math.round(100 * (renderedScale || options.resolution))}% of the window`];
        },
        dispose() {
          controls.enabled = true;
          renderer.debug.onShaderError = previousHandler || null;
          if (target) target.dispose();
          material.dispose();
          blitMaterial.dispose();
        },
      };
    },
  });
})();
