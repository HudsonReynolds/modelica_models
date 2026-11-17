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
    // Compact binary (Examples.BinaryMerger): p = {m1, m2, R1, R2, tMerge, Mf, aF, Rremnant}; the two bodies
    // are the scene's first two tracks. Before tMerge (in ct) the superposed Kerr-Schild metric of the two
    // masses; after it the Kerr remnant (mass Mf, spin length aF about +z), a star of radius Rremnant if > 0.
    binary: (p) => ({ rs: 0, a: 0, q: 0, binary: { m1: p[0], m2: p[1], R1: p[2] || 0, R2: p[3] || 0, tMerge: p[4] ?? 1e30,
      Mf: p[5] ?? (p[0] + p[1]), aF: p[6] || 0, Rrem: p[7] || 0 } }),
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
    // Partial re-tracing with a still camera (see refine()):
    //   uMode 0: everything, bodies included (previews while the camera moves);
    //   uMode 1: the static scene without bodies;
    //   uMode 3: the body mask: alpha 2 where the ray passes through the tube around the stretch of path the
    //            bodies are on (uTube, a 3D occupancy grid over the box uTubeMin + [0, uTubeSize]), else 1;
    //            only those pixels can show a body while it stays on that stretch;
    //   uMode 2: only pixels marked in uMask (3 x 3 dilated) are traced, with bodies; outputs (colour, 1) where
    //            the ray hits a body and 0 elsewhere, to be laid over the static image.
    uniform int uMode;
    uniform highp sampler3D uTube;
    uniform vec3 uTubeMin, uTubeSize;
    uniform float uBodyReach;
    uniform sampler2D uMask;
    // Binary (two bodies before a merger): superposed Kerr-Schild metric of two non-spinning masses uM1, uM2
    // at uX1, uX2 (the bodies' positions now); uR1, uR2 > 0 are star radii (neutron stars), 0 black holes.
    // After the merger uBinary = 0 and the remnant is the Kerr metric of uRs, uA, a star of radius uRrem
    // if it did not collapse.
    uniform int uBinary;
    uniform vec3 uX1, uX2, uStarColour1, uStarColour2, uStarColourRem;
    uniform float uM1, uM2, uR1, uR2, uRrem;
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
    // The binary: g = eta + f1 k1 k1 + f2 k2 k2 with f_a = 2 m_a / r_a, k_a = (1, n_a). The sum has no
    // closed-form inverse, so u = g^{-1} p is a 4 x 4 inverse; then dx/d lambda = u and
    // dp_i/d lambda = 1/2 u^mu u^nu d_i g_mu nu, with d_i f_a = -2 m_a n_i / r_a^2 and
    // d_i k_a,j = (delta_ij - n_i n_j) / r_a.
    mat4 binaryMetric(vec3 x, out vec4 k1, out vec4 k2, out float f1, out float f2, out float r1, out float r2) {
      vec3 d1 = x - uX1, d2 = x - uX2;
      r1 = max(length(d1), 1e-5); r2 = max(length(d2), 1e-5);
      k1 = vec4(1.0, d1 / r1); k2 = vec4(1.0, d2 / r2);
      f1 = 2.0 * uM1 / r1; f2 = 2.0 * uM2 / r2;
      return mat4(-1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0)
        + f1 * outerProduct(k1, k1) + f2 * outerProduct(k2, k2);
    }
    void derivsBinary(vec3 x, vec4 p, out vec3 vel, out vec3 frc) {
      vec4 k1, k2; float f1, f2, r1, r2;
      mat4 G = binaryMetric(x, k1, k2, f1, f2, r1, r2);
      vec4 u = inverse(G) * p;
      vel = u.yzw;
      float a1 = dot(k1, u), a2 = dot(k2, u);
      frc = -uM1 / (r1 * r1) * a1 * a1 * k1.yzw + f1 * a1 * (u.yzw - k1.yzw * dot(k1.yzw, u.yzw)) / r1
            - uM2 / (r2 * r2) * a2 * a2 * k2.yzw + f2 * a2 * (u.yzw - k2.yzw * dot(k2.yzw, u.yzw)) / r2;
    }
    // The metric at x (either case), for the photon's initial momentum.
    mat4 metricAt(vec3 x) {
      if (uBinary == 1) { vec4 k1, k2; float f1, f2, r1, r2; return binaryMetric(x, k1, k2, f1, f2, r1, r2); }
      float f; vec4 k; ks(x, f, k);
      return mat4(-1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0) + f * outerProduct(k, k);
    }
    // Hamilton's equations for H = 1/2 g^{mu nu} p_mu p_nu, g^{-1} = eta - f K K, K^mu = (-k_0, k_i):
    // vel = dx/d lambda = dH/dp, frc = dp/d lambda = -dH/dx.
    void derivs(vec3 x, vec4 p, out vec3 vel, out vec3 frc) {
      if (uBinary == 1) { derivsBinary(x, p, vel, frc); return; }
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

    bool inTube(vec3 x) {
      vec3 tc = (x - uTubeMin) / uTubeSize;
      return all(greaterThanEqual(tc, vec3(0.0))) && all(lessThanEqual(tc, vec3(1.0))) && texture(uTube, tc).r > 0.5;
    }
    void main() {
      if (uMode == 2) {
        ivec2 ms = textureSize(uMask, 0);
        ivec2 mc = ivec2(vUv * vec2(ms));
        float m = 0.0;
        for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++)
          m = max(m, texelFetch(uMask, clamp(mc + ivec2(dx, dy), ivec2(0), ms - 1), 0).a);
        // Written as zero rather than discarded, so the first sample can overwrite the whole layer.
        if (m < 1.5) { gl_FragColor = vec4(0.0); return; }
      }
      vec2 uv = vUv * 2.0 - 1.0 + uJitter;
      vec3 d = normalize(uCamFwd + uv.x * uTanHalfFov * uAspect * uCamRight + uv.y * uTanHalfFov * uCamUp);
      // Null momentum at the camera: contravariant K = (P0, d) with g(K, K) = 0, past-directed so the ray
      // runs back along the light that reached the camera; then lowered, p = g K.
      mat4 Gc = metricAt(uCamPos);
      float A = Gc[0][0], B = 2.0 * dot(vec3(Gc[1][0], Gc[2][0], Gc[3][0]), d);
      float C = dot(d, mat3(Gc[1].yzw, Gc[2].yzw, Gc[3].yzw) * d);
      float P0 = (-B + sqrt(max(B * B - 4.0 * A * C, 0.0))) / (2.0 * A);
      vec4 p = Gc * vec4(P0, d);
      float puObs = p.x / sqrt(max(-A, 1e-6));   // static observer at the camera

      vec3 x = uCamPos;
      vec3 colour = vec3(0.0);
      float transmit = 1.0;
      bool done = false, hit = false;
      bool tube = uMode == 3 && inTube(x);
      vec3 v = velocity(x, p);
      for (int i = 0; i < 4096; i++) {
        if (i >= uMaxSteps) break;
        float r = blR(x);
        // Step size from the distance to the nearest horizon or surface, with a floor proportional to the
        // distance from that body's centre (rays winding onto a horizon still make progress).
        float rNear = r - uHorizon, rCentre = r;
        if (uBinary == 1) {
          // Each body: a horizon at r_a = 2 m_a (as for one Kerr-Schild hole) or a star surface.
          float r1 = length(x - uX1), r2 = length(x - uX2);
          float h1 = uR1 > 0.0 ? uR1 : 2.0 * uM1 * 1.01, h2 = uR2 > 0.0 ? uR2 : 2.0 * uM2 * 1.01;
          if (r1 < h1 || r2 < h2) {
            bool first = r1 - h1 < r2 - h2;
            if ((first ? uR1 : uR2) > 0.0) {
              vec3 n = normalize(x - (first ? uX1 : uX2));
              colour += transmit * (first ? uStarColour1 : uStarColour2) * (0.35 + 0.65 * max(dot(n, -normalize(v)), 0.0));
            }
            done = true; hit = true; break;
          }
          rNear = min(r1 - h1, r2 - h2);
          rCentre = r1 - h1 < r2 - h2 ? r1 : r2;
        } else if (uRrem > 0.0 && r < uRrem) {
          vec3 n = normalize(x);
          colour += transmit * uStarColourRem * (0.35 + 0.65 * max(dot(n, -normalize(v)), 0.0));
          done = true; hit = true; break;
        } else if (uRrem > 0.0) rNear = r - uRrem;
        if (uHorizon > 0.0 && r < uHorizon * 1.01) { done = true; break; }
        if (r > uEscape && dot(x, v) > 0.0) { colour += transmit * sky(normalize(v)); done = true; break; }
        // Outside every body's reach and the photon orbits, an outgoing ray never comes back to a body.
        if (uMode >= 2 && r > uBodyReach && dot(x, v) > 0.0) { done = true; break; }
        float ds = uStepFrac * max(rNear, 0.02 * rCentre + 1e-4);
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
        if (uMode == 3 && inTube(xn)) { tube = true; break; }
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
          if (uMode != 3) colour += transmit * diskColour(mix(x, xn, t), mix(p, pn, t), puObs);
          transmit *= 0.12;
        }
        if (hitBody >= 0) {
          vec3 hp = x + sd * hitT;
          vec3 n = normalize(hp - uBody[hitBody].xyz);
          float lit = 0.35 + 0.65 * max(dot(n, -sd), 0.0);
          colour += transmit * uBodyColor[hitBody] * lit * 1.4;
          done = true; hit = true;
          break;
        }
        x = xn; p = pn;
        v = velocity(x, p);
        if (transmit < 0.02) { done = true; break; }
      }
      // Out of steps: a ray still near the hole is winding onto the horizon (traced backward in time it
      // approaches the past horizon only asymptotically in ingoing Kerr-Schild coordinates), so it is
      // captured; a ray still far out is given the sky in its current direction.
      bool nearHole = uBinary == 1
        ? (uR1 == 0.0 && length(x - uX1) < 8.0 * uM1) || (uR2 == 0.0 && length(x - uX2) < 8.0 * uM2)
        : blR(x) <= 4.0 * max(uHorizon, 1e-6);
      if (!done && !nearHole) colour += transmit * sky(normalize(v));
      // Linear radiance with alpha 1: samples are summed in a float buffer and averaged by alpha at display.
      if (uMode == 2) gl_FragColor = hit ? vec4(colour, 1.0) : vec4(0.0);
      else gl_FragColor = vec4(colour, tube ? 2.0 : 1.0);
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
      // No disk by default without a mass to orbit (Minkowski scenes such as the twin paradox).
      { id: 'disk', label: 'accretion disk', type: 'checkbox', value: (data) => (data && data.metric && FAMILIES[data.metric.name] && FAMILIES[data.metric.name](data.metric.p || []).rs > 0 ? 1 : 0) },
      { id: 'diskOuter', label: 'disk outer (rs)', type: 'range', min: 3, max: 30, step: 0.5, value: 8 },
      { id: 'shift', label: 'Doppler + redshift', type: 'checkbox', value: 1 },
      { id: 'sky', label: 'sky', type: 'select', choices: [{ value: 0, label: 'stars' }, { value: 1, label: 'stars + grid' }, { value: 2, label: 'checker' }], value: 1 },
      { id: 'bodySize', label: 'body size (rs)', type: 'range', min: 0.05, max: 1.5, step: 0.05, value: 0.25 },
      { id: 'camera', label: 'camera', type: 'select', value: 0,
        choices: (data) => [{ value: 0, label: 'orbit (drag)' }].concat(data.trackSpecs.map((t, i) => ({ value: i + 1, label: 'ride with ' + (t.name || 'body ' + (i + 1)) }))) },
      { id: 'paths', label: 'overlay paths', type: 'checkbox', value: 0, help: 'Draw the trajectories on top, without lensing' },
      { id: 'exposure', label: 'exposure', type: 'range', min: 0.2, max: 4, step: 0.1, value: 1.0 },
    ],
    available(data) {
      if (!data.metric) return 'needs a metric record (Visualization.writeMetric, or SceneRecorder metricName)';
      if (!FAMILIES[data.metric.name]) return `the ray tracer does not know the metric "${data.metric.name}"`;
      if (data.metric.name === 'binary' && data.trackSpecs.length < 2) return 'a binary metric needs the two bodies as the first two tracks';
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
      const cLight = data.header.c || 1;
      // A binary: the remnant's horizon and ISCO (it becomes a Kerr hole at tMerge).
      const bin = params.binary || null;
      let remnantHorizon = 0, remnantIsco = 0;
      if (bin) {
        const Mf = bin.Mf, chi = Math.min(Math.abs(bin.aF) / Mf, 0.999);
        remnantHorizon = bin.Rrem > 0 ? 0 : Mf + Math.sqrt(Math.max(Mf * Mf - bin.aF * bin.aF, 0));
        const z1 = 1 + Math.cbrt(1 - chi * chi) * (Math.cbrt(1 + chi) + Math.cbrt(1 - chi));
        const z2 = Math.sqrt(3 * chi * chi + z1 * z1);
        remnantIsco = Mf * (3 + z2 - Math.sqrt((3 - z1) * (3 + z1 + 2 * z2)));
      }
      const kind = (R) => (R > 0 ? `a star of radius ${ctx.fmt(R)}` : 'a black hole');
      const notes = bin
        ? [`Binary: ${bin.R1 > 0 && bin.R2 > 0 ? 'two stars' : bin.R1 > 0 || bin.R2 > 0 ? `${kind(bin.R1)} and ${kind(bin.R2)}` : 'two black holes'} of masses ${ctx.fmt(bin.m1)} and ${ctx.fmt(bin.m2)}${bin.R1 > 0 && bin.R2 > 0 ? ` (radii ${ctx.fmt(bin.R1)}, ${ctx.fmt(bin.R2)})` : ''}, merging at ct = ${bin.tMerge < 1e29 ? ctx.fmt(bin.tMerge) : '(not reached)'} into ${bin.Rrem > 0 ? kind(bin.Rrem) : 'a Kerr black hole'} of mass ${ctx.fmt(bin.Mf)}, a/M = ${ctx.fmt(bin.aF / bin.Mf)}.`,
          'Before the merger the metric is the superposed Kerr-Schild metric of the two masses at their positions now (an approximation, not a solution of the field equations); each frame is a snapshot (no light-travel delay or motion of the bodies while light crosses the scene). After it, the remnant\'s Kerr metric.']
        : [`${data.metric.name}: rs = ${ctx.fmt(params.rs)}, a = ${ctx.fmt(params.a)}, rQ = ${ctx.fmt(params.q)}; horizon r+ = ${ctx.fmt(horizon)}, disk from the ISCO r = ${ctx.fmt(isco)}.`,
          'Snapshot of a stationary metric: bodies are drawn where they are now (no light-travel delay); the camera is a static observer.'];
      if (params.a !== 0) notes.push('Boyer-Lindquist positions are mapped to Kerr-Schild x + iy = (r + ia) sin(theta) e^(i phi), ignoring the small azimuth shift.');
      ctx.setNotes(notes.join('\n'));
      let binaryState = '';

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
        uMode: { value: 0 }, uTube: { value: null }, uTubeMin: { value: new THREE.Vector3() },
        uTubeSize: { value: new THREE.Vector3(1, 1, 1) }, uBodyReach: { value: 1e30 }, uMask: { value: null },
        uBinary: { value: 0 }, uX1: { value: new THREE.Vector3() }, uX2: { value: new THREE.Vector3() },
        uM1: { value: 0 }, uM2: { value: 0 }, uR1: { value: 0 }, uR2: { value: 0 }, uRrem: { value: 0 },
        uStarColour1: { value: new THREE.Color() }, uStarColour2: { value: new THREE.Color() }, uStarColourRem: { value: new THREE.Color(1.1, 1.2, 1.4) },
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
          uniform sampler2D tAccum, tPreview, tDyn; uniform float uUseAccum, uExposure, uDynN;
          void main() {
            vec4 a = texture2D(tAccum, vUv);
            vec4 pv = texture2D(tPreview, vUv);
            // (The preview's alpha carries the body mask, not a sample count.)
            vec3 c = (uUseAccum > 0.5 && a.a > 0.5) ? a.rgb / a.a : pv.rgb;
            // Bodies traced separately: the mean colour of the samples that hit a body, weighted by the
            // fraction of samples that did (anti-aliased edges).
            if (uDynN > 0.5) {
              vec4 d = texture2D(tDyn, vUv);
              if (d.a > 1e-4) c = mix(c, d.rgb / d.a, clamp(d.a / uDynN, 0.0, 1.0));
            }
            c = vec3(1.0) - exp(-c * uExposure);
            gl_FragColor = vec4(pow(c, vec3(1.0 / 2.2)), 1.0);
          }`,
        uniforms: { tAccum: { value: null }, tPreview: { value: null }, tDyn: { value: null }, uUseAccum: { value: 0 }, uExposure: { value: 1 }, uDynN: { value: 0 } },
        depthTest: false, depthWrite: false,
      });
      const blitScene = new THREE.Scene();
      blitScene.add(new THREE.Mesh(new THREE.PlaneGeometry(2, 2), displayMaterial));
      // Half-float targets hold linear radiance and sums of samples; without float rendering only previews are used.
      const floatOk = renderer.capabilities.isWebGL2 && renderer.extensions.has('EXT_color_buffer_float');
      const texType = floatOk ? THREE.HalfFloatType : THREE.UnsignedByteType;
      const makeTarget = (w, h) => new THREE.WebGLRenderTarget(w, h, { type: texType, minFilter: THREE.LinearFilter, magFilter: THREE.LinearFilter, depthBuffer: false });
      let preview = null, accum = null, dyn = null;
      // Refinement state.
      let lastKey = '', lastChange = 0, level = 0, sample = 0, band = 0, bands = 1;
      let pixelBudget = 350000, lastFrame = performance.now(), status = '';
      // Body layer state: the mask is valid after the still-camera preview; the layer restarts when a body moves.
      let maskReady = false, lastBodyKey = '', bodyChange = 0, dynSample = 0, dynFull = false, dynScale = 0.5, dynStatus = '';
      let maskRT = null, lastMask = 0, nMasks = 0;
      // Path length (in units of rs) the tube covers ahead of each body; adapted so the mask is redrawn every
      // second or so, keeping it small without redrawing it every frame.
      let ahead = 4;
      let bodyCount = 0, shownIndex = [];
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
          const size = bin ? Math.max([bin.R1, bin.R2][options.camera - 1] || 0, 2 * [bin.m1, bin.m2][options.camera - 1]) : options.bodySize * unit;
          const out = (P.length() > 1e-9 ? P.clone().normalize() : new THREE.Vector3(1, 0, 0)).multiplyScalar(3 * size + (bin ? 2 * remnantHorizon : 0));
          camera.position.copy(physToWorld(P.clone().add(out)));
          camera.up.set(0, 1, 0);
          camera.lookAt(0, 0, 0);
          controls.enabled = false;
        } else {
          controls.enabled = !ctx.rendering;   // an offline render keeps the camera fixed
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
        // A binary: the first two tracks are the two masses (always traced; they are the spacetime), then the
        // remnant from the merger on. Its state is part of the camera key: everything is re-traced when it moves.
        let binKey = '';
        if (bin) {
          const merged = cLight * ctx.clock >= bin.tMerge;
          uniforms.uBinary.value = merged ? 0 : 1;
          uniforms.uM1.value = bin.m1; uniforms.uM2.value = bin.m2;
          uniforms.uR1.value = bin.R1; uniforms.uR2.value = bin.R2;
          uniforms.uX1.value.copy(samplePhys(bodies[0].sample.x));
          uniforms.uX2.value.copy(samplePhys(bodies[1].sample.x));
          const star = new THREE.Color(0.85, 0.92, 1.0);
          uniforms.uStarColour1.value.copy(star).lerp(bodies[0].color, 0.3).multiplyScalar(1.5);
          uniforms.uStarColour2.value.copy(star).lerp(bodies[1].color, 0.3).multiplyScalar(1.5);
          uniforms.uRs.value = merged ? 2 * bin.Mf : 0;
          uniforms.uA.value = merged ? bin.aF : 0;
          uniforms.uHorizon.value = merged ? remnantHorizon : 0;
          uniforms.uRrem.value = merged ? bin.Rrem : 0;
          uniforms.uDiskIn.value = merged ? remnantIsco : unit;
          const sep = uniforms.uX1.value.distanceTo(uniforms.uX2.value);
          binaryState = merged ? `merged: ${bin.Rrem > 0 ? 'star' : 'Kerr hole'} M = ${ctx.fmt(bin.Mf)}, a/M = ${ctx.fmt(bin.aF / bin.Mf)}` : `inspiral: separation ${ctx.fmt(sep)}`;
          binKey = [merged, uniforms.uX1.value.toArray().map(v => v.toPrecision(6)).join(','), uniforms.uX2.value.toArray().map(v => v.toPrecision(6)).join(',')].join('/');
        }
        // Bodies hidden in the legend are not traced (a binary's bodies are the spacetime itself).
        shownIndex = bin ? [] : bodies.map((b, i) => (b.hidden ? -1 : i)).filter(i => i >= 0).slice(0, 4);
        const count = shownIndex.length;
        bodyCount = count;
        uniforms.uBodyCount.value = count;
        for (let i = 0; i < count; i++) {
          const P = samplePhys(bodies[shownIndex[i]].sample.x);
          uniforms.uBody.value[i].set(P.x, P.y, P.z, options.bodySize * unit);
          uniforms.uBodyColor.value[i].copy(bodies[shownIndex[i]].color);
        }
        // What changes the image: the camera key (camera, options, window, which bodies are shown) restarts
        // everything; the body key (where they are) only the body layer. The camera is rounded so that motion
        // far below a pixel (the tail of an easing camera) does not restart refinement.
        const { exposure, paths, ...traced } = options;
        return {
          cam: [e.map(v => v.toPrecision(5)).join(','), shownIndex.join(','), JSON.stringify(traced), renderer.domElement.width, renderer.domElement.height, binKey].join('|'),
          bodies: uniforms.uBody.value.slice(0, count).map(b => b.toArray().join(',')).join(';'),
        };
      }

      // The tube: every place a shown body can be, over the whole recording, widened by the body radius, half
      // a step of the integrator there (rays are tested only at step ends) and a margin for the slightly
      // different paths of the finer refinement steps, marked in a 3D occupancy grid. A pixel whose ray never
      // enters it shows no body at any time, so it never needs re-tracing while the camera holds still.
      let tubeTex = null, tubeInfo = '', tubeWindow = null;
      // Sample index ranges [from, to] per shown track: from just behind the body to `ahead` rs in front.
      function tubeWindowNow(bodies) {
        return shownIndex.map((ti) => {
          const smp = data.trackSpecs[ti].samples, cur = bodies[ti].sample.index || 0;
          let from = cur, to = cur, len = 0;
          while (from > 0 && len < 0.5 * unit) { len += samplePhys(smp[from].x).distanceTo(samplePhys(smp[from - 1].x)); from--; }
          len = 0;
          while (to < smp.length - 1 && len < ahead * unit) { len += samplePhys(smp[to].x).distanceTo(samplePhys(smp[to + 1].x)); to++; }
          return { ti, from, to: Math.min(to + 1, smp.length - 1), cur };
        });
      }
      // The window still holds every body (with a sample to spare ahead, unless at the end of its track).
      function windowHolds(win, bodies) {
        return win && win.length === shownIndex.length && win.every((w, k) => {
          const cur = bodies[shownIndex[k]].sample.index || 0, last = data.trackSpecs[w.ti].samples.length - 1;
          return w.ti === shownIndex[k] && cur >= w.from && (cur < w.to - 1 || w.to === last);
        });
      }
      function buildTube(win) {
        tubeWindow = win;
        const R = options.bodySize * unit;
        const widen = (P) => { const rr = P.length(); return R + 0.5 * options.accuracy * Math.max(rr - horizon, 0.02 * rr) + 0.25 * unit; };
        // Positions along each shown track, interpolated in the recorded coordinates (as the bodies move).
        const pts = [];
        for (const w of win) {
          const smp = data.trackSpecs[w.ti].samples;
          for (let k = w.from; k <= w.to; k++) {
            const a = smp[k], b = smp[Math.min(k + 1, w.to)];
            const Pa = samplePhys(a.x), Pb = samplePhys(b.x);
            const n = k < w.to ? Math.max(1, Math.ceil(Pa.distanceTo(Pb) / (0.25 * R))) : 1;
            for (let j = 0; j < n; j++) {
              const f = j / n;
              pts.push(samplePhys(a.x.map((v, m) => v + f * (b.x[m] - v))));
            }
          }
        }
        if (!pts.length) return false;
        const lo = new THREE.Vector3(Infinity, Infinity, Infinity), hi = new THREE.Vector3(-Infinity, -Infinity, -Infinity);
        let reach = 0;
        for (const P of pts) {
          const w = widen(P);
          lo.min(P.clone().subScalar(w)); hi.max(P.clone().addScalar(w));
          reach = Math.max(reach, P.length() + w);
        }
        const size = hi.clone().sub(lo);
        const vox = Math.max(size.x, size.y, size.z) / 96;
        const nx = Math.max(2, Math.ceil(size.x / vox)), ny = Math.max(2, Math.ceil(size.y / vox)), nz = Math.max(2, Math.ceil(size.z / vox));
        const cell = new THREE.Vector3(size.x / nx, size.y / ny, size.z / nz);
        const half = 0.5 * cell.length();
        const grid = new Uint8Array(nx * ny * nz);
        for (const P of pts) {
          const rad = widen(P) + half;
          const i0 = Math.max(0, Math.floor((P.x - rad - lo.x) / cell.x)), i1 = Math.min(nx - 1, Math.floor((P.x + rad - lo.x) / cell.x));
          const j0 = Math.max(0, Math.floor((P.y - rad - lo.y) / cell.y)), j1 = Math.min(ny - 1, Math.floor((P.y + rad - lo.y) / cell.y));
          const k0 = Math.max(0, Math.floor((P.z - rad - lo.z) / cell.z)), k1 = Math.min(nz - 1, Math.floor((P.z + rad - lo.z) / cell.z));
          for (let k = k0; k <= k1; k++) for (let j = j0; j <= j1; j++) for (let i = i0; i <= i1; i++) {
            const dx = lo.x + (i + 0.5) * cell.x - P.x, dy = lo.y + (j + 0.5) * cell.y - P.y, dz = lo.z + (k + 0.5) * cell.z - P.z;
            if (dx * dx + dy * dy + dz * dz <= rad * rad) grid[(k * ny + j) * nx + i] = 255;
          }
        }
        if (tubeTex) tubeTex.dispose();
        tubeTex = new THREE.Data3DTexture(grid, nx, ny, nz);
        tubeTex.format = THREE.RedFormat;
        tubeTex.type = THREE.UnsignedByteType;
        tubeTex.minFilter = tubeTex.magFilter = THREE.NearestFilter;
        tubeTex.unpackAlignment = 1;
        tubeTex.needsUpdate = true;
        uniforms.uTube.value = tubeTex;
        uniforms.uTubeMin.value.copy(lo);
        uniforms.uTubeSize.value.copy(size);
        // Photon orbits reach r = 2 rs (retrograde, extreme Kerr); beyond that and the tube, outgoing rays escape.
        uniforms.uBodyReach.value = Math.max(reach, 2 * params.rs);
        tubeInfo = `${nx} x ${ny} x ${nz}`;
        return true;
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
      function trace(r, targetRT, mat, stepFrac, maxSteps, jx, jy, y0, h, mode = 0) {
        uniforms.uMode.value = mode;
        uniforms.uBodyCount.value = mode === 1 || mode === 3 ? 0 : bodyCount;
        uniforms.uStepFrac.value = stepFrac;
        uniforms.uMaxSteps.value = maxSteps;
        uniforms.uJitter.value.set(jx, jy);
        quad.material = mat;
        // three.js re-uploads a material's uniforms only when the material changes between draws; two passes of
        // the same material in a row (a static band, then the body layer) need it forced.
        mat.uniformsNeedUpdate = true;
        if (h !== undefined) { targetRT.scissor.set(0, y0, targetRT.width, h); targetRT.scissorTest = true; }
        else targetRT.scissorTest = false;
        // autoClear would clear the target (or the scissored band) before every pass, so additive samples would
        // never sum; previews overwrite every pixel anyway, and the accumulation buffer is cleared explicitly.
        const autoClear = r.autoClear;
        r.autoClear = false;
        r.setRenderTarget(targetRT);
        r.render(quadScene, quadCamera);
        r.setRenderTarget(null);
        r.autoClear = autoClear;
        targetRT.scissorTest = false;
      }
      // One frame of progressive refinement: a coarse preview while anything changes, a sharper preview once
      // still, then full-resolution samples with finer steps and sub-pixel jitter, traced in horizontal bands
      // that fit a per-frame pixel budget and summed into the accumulation buffer.
      //
      // With a still camera only the bodies change. With float targets the passes from the sharper preview on
      // trace the scene without bodies, and that preview also marks the pixels whose rays cross the bodies'
      // tube. A body layer then re-traces only those pixels when a body moves (one sample at an adaptive
      // resolution while moving; full resolution, accumulating samples, once the bodies stop), laid over the
      // refined static image, which keeps refining while the bodies move.
      let frameTime = 0;
      const slowFrame = () => (ctx.rendering ? 250 : 70);
      function refine(r, keys) {
        const now = performance.now();
        frameTime = now - lastFrame;
        lastFrame = now;
        const key = floatOk ? keys.cam : keys.cam + '#' + keys.bodies;   // no float targets: re-trace everything
        if (key !== lastKey) { lastKey = key; lastChange = now; level = 0; maskReady = false; }
        // An offline render does not wait for stillness: it traces each frame to convergence.
        const still = ctx.rendering || now - lastChange > 150;
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
          trace(r, preview, material, options.accuracy, Math.round(options.steps), 0, 0, undefined, undefined, floatOk ? 1 : 0);
          maskReady = false;
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
          // Adapt the per-frame budget to how long the previous refinement frame took (an offline render
          // accepts slower frames: nobody is interacting).
          if (frameTime > slowFrame()) pixelBudget = Math.max(60000, pixelBudget * 0.7);
          bands = Math.max(1, Math.ceil(w * h / pixelBudget));
        } else if (band === 0) {
          if (frameTime > slowFrame()) pixelBudget = Math.max(60000, pixelBudget * 0.7);
          else if (frameTime < 0.4 * slowFrame()) pixelBudget = Math.min(ctx.rendering ? 1.6e7 : 4e6, pixelBudget * 1.25);
          bands = Math.max(1, Math.ceil(w * h / pixelBudget));
        }
        const bandH = Math.ceil(h / bands), y0 = band * bandH;
        // Halton (2, 3) jitter within the pixel; the first sample is centred.
        const jx = sample === 0 ? 0 : (halton(sample, 2) - 0.5) * 2 / w;
        const jy = sample === 0 ? 0 : (halton(sample, 3) - 0.5) * 2 / h;
        trace(r, accum, addMaterial, options.accuracy * 0.6, Math.round(options.steps * 1.5), jx, jy, y0, Math.min(bandH, h - y0), 1);
        band++;
        if (band >= bands) { band = 0; sample++; }
        if (sample >= options.samples) level = 3;
        status = level === 3 ? `refined: ${sample} samples at ${w} x ${h}` : `refining: sample ${sample + 1} of ${options.samples} at ${w} x ${h}`;
        return true;
      }
      // The body layer (after refine, in the same frame). Returns the number of samples it holds.
      function refineBodies(r, keys, bodies) {
        if (!floatOk || bodyCount === 0 || level < 2) { dynStatus = ''; return 0; }
        const now = performance.now();
        // (Re)draw the mask when the camera settled or a body is about to leave the stretch of path it covers.
        if (!maskReady || !windowHolds(tubeWindow, bodies)) {
          if (maskReady) {
            const dt = now - lastMask;
            if (dt < 700) ahead = Math.min(ahead * 1.6, 400);
            else if (dt > 2500) ahead = Math.max(ahead / 1.3, 2);
          }
          if (!buildTube(tubeWindowNow(bodies))) { dynStatus = ''; return 0; }
          const el = renderer.domElement;
          maskRT = ensure(maskRT, Math.max(16, Math.round(el.clientWidth * 0.5)), Math.max(16, Math.round(el.clientHeight * 0.5)));
          trace(r, maskRT, material, options.accuracy, Math.round(options.steps), 0, 0, undefined, undefined, 3);
          maskReady = true; lastMask = now; nMasks++;
          lastBodyKey = '';   // the layer is redrawn with the new mask
        }
        if (keys.bodies !== lastBodyKey) {
          if (lastBodyKey && dynSample > 0 && !dynFull) {
            // Moving: keep the body layer within a frame budget by adapting its resolution.
            if (frameTime > 45) dynScale = Math.max(0.25, dynScale * 0.85);
            else if (frameTime < 22) dynScale = Math.min(1, dynScale * 1.1);
          }
          // An offline render goes straight to full-resolution samples (the bodies are still in each frame).
          lastBodyKey = keys.bodies; bodyChange = now; dynSample = 0; dynFull = !!ctx.rendering;
        }
        if (!dynFull && dynSample > 0 && now - bodyChange > 150) { dynFull = true; dynSample = 0; }
        const target = dynFull ? Math.max(1, Math.round(options.samples)) : 1;
        if (dynSample < target) {
          const [W, H] = passSize(2), sc = dynFull ? 1 : dynScale;
          const w = Math.max(16, Math.round(W * sc)), h = Math.max(16, Math.round(H * sc));
          dyn = ensure(dyn, w, h);
          uniforms.uMask.value = maskRT.texture;
          const jx = dynSample === 0 ? 0 : (halton(dynSample, 2) - 0.5) * 2 / w;
          const jy = dynSample === 0 ? 0 : (halton(dynSample, 3) - 0.5) * 2 / h;
          // The first sample overwrites the layer (every pixel is written), later ones add to it.
          trace(r, dyn, dynSample === 0 ? material : addMaterial, options.accuracy * (dynFull ? 0.6 : 1), Math.round(options.steps * (dynFull ? 1.5 : 1)), jx, jy, undefined, undefined, 2, 0);
          // Unbound again: a texture bound while its own target is drawn into would be a feedback loop.
          uniforms.uMask.value = null;
          dynSample++;
          dynStatus = `bodies: ${dynFull ? `sample ${dynSample} of ${target}` : 'moving'} at ${w} x ${h}, only near their paths (${ahead.toFixed(0)} rs ahead; mask ${nMasks})`;
        }
        return dynSample;
      }

      return {
        group: new THREE.Group(),
        hidesBodies: true,
        place(x) { return physToWorld(samplePhys(x)); },
        camera() {
          // A binary: far enough to see the whole initial orbit (sceneSize is twice the largest orbital radius).
          const R = bin ? Math.max(2.4 * ctx.sceneSize, 30 * unit) : Math.min(Math.max(14 * unit, 0.9 * ctx.extent), 60 * unit);
          // Slightly above the disk plane, so the disk's far side is lensed over the shadow; a binary from higher
          // up, so the orbit is seen opening out.
          const incl = THREE.MathUtils.degToRad(bin ? 62 : 80);
          return { position: new THREE.Vector3(R * Math.sin(incl), R * Math.cos(incl), R * 0.15), target: new THREE.Vector3(0, 0, 0) };
        },
        render(r, cam, scene3) {
          const keys = syncUniforms();
          const useAccum = refine(r, keys);
          const nDyn = refineBodies(r, keys, ctx.bodiesNow());
          displayMaterial.uniforms.tAccum.value = accum ? accum.texture : preview.texture;
          displayMaterial.uniforms.tPreview.value = preview.texture;
          displayMaterial.uniforms.uUseAccum.value = useAccum && accum ? 1 : 0;
          displayMaterial.uniforms.tDyn.value = nDyn > 0 ? dyn.texture : null;
          displayMaterial.uniforms.uDynN.value = nDyn;
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
        // For offline rendering: the frame is final (all samples, including the body layer's).
        converged() {
          if (level !== 3) return false;
          if (floatOk && bodyCount > 0 && maskReady) return dynFull && dynSample >= Math.max(1, Math.round(options.samples));
          return true;
        },
        readout() {
          return [bin ? binaryState : `horizon r+ ${ctx.fmt(horizon)}, ISCO ${ctx.fmt(isco)}`, status].concat(dynStatus ? [dynStatus] : []);
        },
        dispose() {
          controls.enabled = true;
          controls.enableDamping = previousDamping;
          renderer.debug.onShaderError = previousHandler || null;
          if (preview) preview.dispose();
          if (accum) accum.dispose();
          if (dyn) dyn.dispose();
          if (maskRT) maskRT.dispose();
          if (tubeTex) tubeTex.dispose();
          material.dispose();
          addMaterial.dispose();
          displayMaterial.dispose();
        },
      };
    },
  });
})();
