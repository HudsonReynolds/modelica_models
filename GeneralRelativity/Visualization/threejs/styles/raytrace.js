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

  // Blackbody colours. CIE 1931 colour-matching functions in the multi-lobe Gaussian fit of Wyman, Sloan and
  // Shirley (JCGT 2(2), 2013); Planck's law integrated against them gives XYZ, then linear sRGB (D65).
  const lobe = (l, mu, s1, s2) => { const t = (l - mu) / (l < mu ? s1 : s2); return Math.exp(-0.5 * t * t); };
  const cieX = (l) => 1.056 * lobe(l, 599.8, 37.9, 31.0) + 0.362 * lobe(l, 442.0, 16.0, 26.7) - 0.065 * lobe(l, 501.1, 20.4, 26.2);
  const cieY = (l) => 0.821 * lobe(l, 568.8, 46.9, 40.5) + 0.286 * lobe(l, 530.9, 16.3, 31.1);
  const cieZ = (l) => 1.217 * lobe(l, 437.0, 11.8, 36.0) + 0.681 * lobe(l, 459.0, 26.0, 13.8);
  // [r, g, b per unit luminance Y, log10 Y] of a blackbody at T kelvin (Y in arbitrary but fixed units).
  function blackbodyRGB(T) {
    let X = 0, Y = 0, Z = 0;
    for (let l = 360; l <= 830; l += 5) {
      const um = l * 1e-3;
      const B = 1 / (Math.pow(um, 5) * Math.expm1(14387.77 / (um * T)));   // second radiation constant hc/k in um K
      X += B * cieX(l); Y += B * cieY(l); Z += B * cieZ(l);
    }
    if (!(Y > 0)) return [0, 0, 0, -300];
    return [(3.2406 * X - 1.5372 * Y - 0.4986 * Z) / Y, (-0.9689 * X + 1.8758 * Y + 0.0415 * Z) / Y, (0.0557 * X - 0.2040 * Y + 1.0570 * Z) / Y, Math.log10(Y)];
  }
  // Mean visible luminance (the units of blackbodyRGB) of a disk with peak temperature Tpeak over its area out to
  // 4 r_in: the disk's brightness unit, so that hot and cool disks show at a similar overall brightness (like a
  // camera's automatic exposure) while the variation across the disk stays physical.
  function diskMeanLuminance(Tpeak) {
    let num = 0, den = 0;
    for (let i = 0; i < 48; i++) {
      const r = 1 + 3 * (i + 0.5) / 48, x = 1 / r;
      const temp = Math.pow(x, 0.75) * Math.pow(1 - Math.sqrt(x), 0.25) / 0.488;
      num += Math.pow(10, blackbodyRGB(Tpeak * temp)[3]) * r;
      den += r;
    }
    return Math.max(num / den, 1e-300);
  }
  // The table the shader interpolates: log10 T from 2.5 (316 K) to 8.5 (3e8 K).
  const BB_N = 81, BB_LO = 2.5, BB_STEP = 6 / (BB_N - 1);
  const BB_TABLE = Array.from({ length: BB_N }, (_, i) => blackbodyRGB(Math.pow(10, BB_LO + i * BB_STEP)));
  // Peak temperature (K) of a thin disk (Shakura and Sunyaev, no torque at the inner edge) around a hole of
  // massSolar solar masses accreting at the fraction `eddington` of the Eddington rate, with its inner edge at
  // rInOverM = r_isco / M: T^4 = 3 G M Mdot / (8 pi sigma r^3) (1 - sqrt(r_in / r)), largest at r = 49/36 r_in,
  // with Mdot = eddington L_Edd / (eta c^2) and the efficiency eta = 1 - E_isco = 1 - sqrt(1 - 2 M / (3 r_isco)).
  // About 4.9e6 K for 10 solar masses at a tenth of Eddington, scaling as M^(-1/4).
  function diskPeakTemperature(massSolar, eddington, rInOverM) {
    const G = 6.6743e-11, c = 2.99792458e8, sigmaSB = 5.670374e-8, sigmaT = 6.6524587e-29, mp = 1.67262192e-27, Msun = 1.98847e30;
    const M = massSolar * Msun;
    const eta = 1 - Math.sqrt(1 - 2 / (3 * rInOverM));
    const Mdot = eddington * 4 * Math.PI * G * M * mp * c / sigmaT / (eta * c * c);
    const r = 49 / 36 * rInOverM * G * M / (c * c);
    return Math.pow(3 * G * M * Mdot / (8 * Math.PI * sigmaSB * r * r * r) * (1 - 6 / 7), 0.25);
  }

  // The disk's own options are listed while the disk is on.
  const diskShown = (data, o) => !o || !!o.disk;
  const familyOf = (data) => (data && data.metric && FAMILIES[data.metric.name] ? FAMILIES[data.metric.name](data.metric.p || []) : null);
  const neutronStars = (data) => { const f = familyOf(data); return !!(f && f.binary && f.binary.R1 > 0 && f.binary.R2 > 0); };
  const mergerScene = (data) => { const f = familyOf(data); return !!(f && f.binary && f.binary.tMerge < 1e29); };
  const smooth = (x, a, b) => { const t = Math.min(Math.max((x - a) / (b - a), 0), 1); return t * t * (3 - 2 * t); };

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
    // Disk emission: peak temperature (K), log10 of the disk's mean visible luminance (the brightness unit), strength of the
    // structure, the shear age of the pattern, and up to two pattern layers (weight, time offset, seed) that
    // cross-fade while the disk turns. uBB: blackbody colours, see planck().
    uniform float uDiskTPeak, uDiskNorm, uDiskDetail, uDiskAge;
    uniform vec4 uDiskLayer[2];
    const int BB_N = ${BB_N};
    const float BB_LO = ${BB_LO.toFixed(4)}, BB_STEP = ${BB_STEP.toFixed(6)};
    uniform vec4 uBB[BB_N];
    uniform vec2 uJitter;   // sub-pixel offset in uv units, for accumulated anti-aliasing
    uniform float uFrameSeed;   // changes with the clock: the gas dither does not stand still in a video
    uniform float uPixelAngle;   // angle across one pixel of the pass being traced, at the centre of the view (rad)
    uniform int uBodyCount;
    uniform vec4 uBody[4];
    uniform vec3 uBodyColor[4];
    // Lighting of the bodies: 0 by the accretion disk (a distant star without one), 1 a distant star in the
    // direction uSunDir, 2 from the camera; uBodyGlow: their own faint glow; uShadowR: the radius around the hole
    // that blocks straight light paths (about the photon sphere).
    uniform int uLightMode;
    uniform vec3 uSunDir;
    uniform float uBodyGlow, uShadowR;
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
    // Tidal shapes of the stars (uTides = 1), as one field D, negative inside and about the distance to the
    // surface (shapeField): the Roche potential of point masses uM1, uM2 at uX1, uX2 in the frame rotating at
    // uOmega2^(1/2) about the centre of mass uCM, below a surface level (uLevel1 on star 1's side of the plane
    // through the inner Lagrange point uL1, uLevel2 on the other), inside the bounding spheres uB1, uB2 that hold
    // the surface; blended (uBlend) into a ball blR = uBall about the origin, the merged remnant. uHeat whitens
    // shocked matter. uRemnant = 1: the remnant ball after the Kerr metric has taken over, until it has
    // collapsed onto the horizon.
    uniform int uTides, uRemnant;
    uniform vec3 uCM, uL1;
    uniform float uOmega2, uLevel1, uLevel2, uB1, uB2, uGFloor, uBlend, uBall, uHeat;
    // Gas as a volume (uDebris): 1 the gas of a neutron-star merger, 2 a thick accretion disk. The merger's gas
    // is a set of particles (see the gas simulation in create) deposited into the 3D texture uGasTex over the
    // box uGasMin + [0, uGasSize]: per voxel the density, and the density of gas still young enough to show
    // the colour of star 1 and of star 2; uGasScale turns it into optical depth per unit length. A thick disk
    // is a torus of thickness H = uDebH r (Gaussian in z) from uDebIn to uDebOut, filled out to uFront with
    // weight uFill. Both orbit a mass uDebM with spin length uDebA.
    uniform int uDebris;
    uniform float uDebM, uDebA, uDebIn, uDebOut, uDebH, uFill, uFront;
    uniform highp sampler3D uGasTex;
    uniform vec3 uGasMin, uGasSize;
    uniform float uGasScale;
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
    // While a merger's remnant forms, the two masses carry the weight uBinW and the remnant's Kerr-Schild term
    // (uRs, uA) the weight uKerrW = 1 - uBinW: the pair's rotating lensing gives way to the spinning hole's
    // without a stop in between, and at uKerrW = 1 the metric is exactly the Kerr branch's.
    uniform float uBinW, uKerrW;
    mat4 binaryMetric(vec3 x, out vec4 k1, out vec4 k2, out float f1, out float f2, out float r1, out float r2) {
      vec3 d1 = x - uX1, d2 = x - uX2;
      r1 = max(length(d1), 1e-5); r2 = max(length(d2), 1e-5);
      k1 = vec4(1.0, d1 / r1); k2 = vec4(1.0, d2 / r2);
      f1 = 2.0 * uBinW * uM1 / r1; f2 = 2.0 * uBinW * uM2 / r2;
      mat4 G = mat4(-1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0)
        + f1 * outerProduct(k1, k1) + f2 * outerProduct(k2, k2);
      if (uKerrW > 0.0) { float f; vec4 k; ks(x, f, k); G += uKerrW * f * outerProduct(k, k); }
      return G;
    }
    void derivsBinary(vec3 x, vec4 p, out vec3 vel, out vec3 frc) {
      vec4 k1, k2; float f1, f2, r1, r2;
      mat4 G = binaryMetric(x, k1, k2, f1, f2, r1, r2);
      vec4 u = inverse(G) * p;
      vel = u.yzw;
      float a1 = dot(k1, u), a2 = dot(k2, u), m1 = uBinW * uM1, m2 = uBinW * uM2;
      frc = -m1 / (r1 * r1) * a1 * a1 * k1.yzw + f1 * a1 * (u.yzw - k1.yzw * dot(k1.yzw, u.yzw)) / r1
            - m2 / (r2 * r2) * a2 * a2 * k2.yzw + f2 * a2 * (u.yzw - k2.yzw * dot(k2.yzw, u.yzw)) / r2;
      if (uKerrW > 0.0) {
        float f; vec4 k; vec3 gf, gk1, gk2, gk3;
        ksGrad(x, f, k, gf, gk1, gk2, gk3);
        float kv = dot(k, u);
        frc += uKerrW * (0.5 * gf * kv * kv + f * kv * (gk1 * u.y + gk2 * u.z + gk3 * u.w));
      }
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
    // The blackbody table at T kelvin (colour per unit luminance, log10 luminance), interpolated in log10 T.
    vec4 bbEntry(float T) {
      float u = clamp((log(max(T, 1.0)) * 0.4342945 - BB_LO) / BB_STEP, 0.0, float(BB_N - 1) - 1e-3);
      int i = int(u);
      return mix(uBB[i], uBB[i + 1], u - float(i));
    }
    // Cube-map faces 0..5 (+x, -x, +y, -y, +z, -z) with coordinates uv in [-1, 1].
    void cubeFace(vec3 d, out int f, out vec2 uv) {
      vec3 a = abs(d);
      if (a.x >= a.y && a.x >= a.z) { f = d.x > 0.0 ? 0 : 1; uv = d.yz / a.x; }
      else if (a.y >= a.z) { f = d.y > 0.0 ? 2 : 3; uv = d.xz / a.y; }
      else { f = d.z > 0.0 ? 4 : 5; uv = d.xy / a.z; }
    }
    vec3 cubeDir(int f, vec2 uv) {
      float s = (f == 1 || f == 3 || f == 5) ? -1.0 : 1.0;
      return normalize(f < 2 ? vec3(s, uv) : f < 4 ? vec3(uv.x, s, uv.y) : vec3(uv, s));
    }
    // Stars are point sources: their images keep the size of a pixel however lensing stretches the sky around
    // them, and lensing changes only their brightness, by the magnification mu (as in James et al., Class. Quantum
    // Grav. 32, 065001 (2015)). sx and sy are the sky-direction steps to the next pixel in x and y, so a sky offset
    // dd is the pixel offset p with [sx sy] p = dd (least squares, through G = [sx sy]^T [sx sy] = (ga gb; gb gc)).
    // At most one star per cell of an n x n grid on each cube face, kept inside its cell (so never cut off), drawn
    // as a Gaussian of 0.55 pixels. Fluxes follow N(>F) ~ F^-1.5 (stars spread evenly in space) from fMin up to
    // 40 fMin; colours are blackbodies of 3000 to 30000 K; galactic raises the density towards the Milky Way band.
    vec3 stars(vec3 dir, float n, float density, float fMin, float galactic, vec3 sx, vec3 sy, vec3 G, float mu) {
      int f; vec2 uv; cubeFace(dir, f, uv);
      vec2 c = floor((uv * 0.5 + 0.5) * n);
      vec3 key = vec3(c, float(f) + 0.37 * n);
      float b = dot(cubeDir(f, (c + 0.5) / n * 2.0 - 1.0), normalize(vec3(0.25, 0.55, 0.8)));
      float dens = density * (1.0 + galactic * (2.5 * exp(-12.0 * b * b) - 0.6));
      float h = hash13(key);
      if (h < 1.0 - dens) return vec3(0.0);
      vec3 dd = cubeDir(f, (c + 0.25 + 0.5 * vec2(hash13(key + 1.7), hash13(key + 3.1))) / n * 2.0 - 1.0) - dir;
      vec2 q = vec2(dot(sx, dd), dot(sy, dd));
      vec2 p = vec2(G.z * q.x - G.y * q.y, G.x * q.y - G.y * q.x) / (G.x * G.z - G.y * G.y);
      float g = exp(-dot(p, p) / (2.0 * 0.55 * 0.55));
      if (g < 1e-3) return vec3(0.0);
      float flux = fMin * pow(max((h - (1.0 - dens)) / dens, 0.004), -0.6667);
      vec3 tint = max(bbEntry(3000.0 * pow(10.0, hash13(key + 5.3))).rgb, 0.0);
      return tint * (flux * mu * g);
    }
    // The sky in direction dir. sx, sy: its change to the next pixel in x and y; camArea: the solid angle of the
    // pixel at the camera (|dFdx(d) x dFdy(d)|); valid false where the neighbours' rays did not reach the sky
    // (then an unlensed pixel of uPixelAngle is assumed).
    vec3 sky(vec3 dir, vec3 sx, vec3 sy, float camArea, bool valid) {
      float px2 = uPixelAngle * uPixelAngle;
      vec3 G = vec3(dot(sx, sx), dot(sx, sy), dot(sy, sy));
      float area = sqrt(max(G.x * G.z - G.y * G.y, 0.0));
      if (!valid || !(area > 1e-4 * px2)) {
        sx = normalize(cross(dir, abs(dir.z) < 0.9 ? vec3(0.0, 0.0, 1.0) : vec3(1.0, 0.0, 0.0))) * uPixelAngle;
        sy = cross(dir, sx);
        G = vec3(px2, 0.0, px2);
        area = camArea = px2;
      }
      float mu = clamp(camArea / area, 0.0, 20.0);   // magnification: solid angle at the camera over that on the sky
      float th = acos(clamp(dir.z, -1.0, 1.0));
      float ph = atan(dir.y, dir.x);
      vec3 col = vec3(0.0006, 0.0008, 0.0016);
      float band = dot(dir, normalize(vec3(0.25, 0.55, 0.8)));
      col += vec3(0.010, 0.009, 0.014) * exp(-10.0 * band * band);
      col += stars(dir, 160.0, 0.10, 0.15, 0.0, sx, sy, G, mu) + stars(dir, 420.0, 0.06, 0.05, 1.0, sx, sy, G, mu);
      if (uSky == 1) {
        // Lines of latitude and longitude every 15 degrees, about 1.5 pixels wide wherever lensing puts them: the
        // angular distance to a line over its change per pixel across it (along e_theta or e_phi).
        float spacing = PI / 12.0;
        float lat = abs(fract(th / spacing + 0.5) - 0.5) * spacing;
        float lon = abs(fract(ph / spacing + 0.5) - 0.5) * spacing * sin(th);
        vec3 eTh = vec3(cos(th) * cos(ph), cos(th) * sin(ph), -sin(th)), ePh = vec3(-sin(ph), cos(ph), 0.0);
        float pxLat = lat / max(length(vec2(dot(eTh, sx), dot(eTh, sy))), 1e-9);
        float pxLon = lon / max(length(vec2(dot(ePh, sx), dot(ePh, sy))), 1e-9);
        col += vec3(0.045, 0.10, 0.17) * (1.0 - smoothstep(0.3, 1.0, min(pxLat, pxLon)));
      } else if (uSky == 2) {
        float c = mod(floor(th / (PI / 12.0)) + floor((ph + PI) / (PI / 12.0)), 2.0);
        col = mix(vec3(0.01, 0.012, 0.025), vec3(0.08, 0.09, 0.13), c) + col * 0.5;
      }
      return col;
    }
    // Visible light of a blackbody at T kelvin as linear sRGB, in units of the disk's mean (uDiskNorm). uBB holds,
    // at log10 T = BB_LO + i BB_STEP, the colour per unit luminance and log10 of the luminance, from Planck's
    // law and the CIE 1931 matching functions (blackbodyRGB in JS).
    vec3 planck(float T) {
      if (log(max(T, 1.0)) * 0.4342945 <= BB_LO) return vec3(0.0);
      vec4 e = bbEntry(T);
      return max(e.rgb, 0.0) * exp2((e.a - uDiskNorm) * 3.3219281);
    }
    float vnoise(vec3 p) {
      vec3 i = floor(p), f = fract(p);
      f = f * f * (3.0 - 2.0 * f);
      return mix(mix(mix(hash13(i), hash13(i + vec3(1.0, 0.0, 0.0)), f.x), mix(hash13(i + vec3(0.0, 1.0, 0.0)), hash13(i + vec3(1.0, 1.0, 0.0)), f.x), f.y),
                 mix(mix(hash13(i + vec3(0.0, 0.0, 1.0)), hash13(i + vec3(1.0, 0.0, 1.0)), f.x), mix(hash13(i + vec3(0.0, 1.0, 1.0)), hash13(i + vec3(1.0, 1.0, 1.0)), f.x), f.y), f.z);
    }
    // Structure of the disk at radius r and azimuth phi, where it turns at omega: streaks along the orbits, sheared
    // into trailing spirals (the pattern is advected for a time uDiskAge plus the layer's offset), rings and gaps
    // in radius, and warm and cool patches. Returns (coverage: the fraction of the surface that glows and blocks
    // light, temperature factor). The angle enters through (cos, sin), so the pattern is periodic in phi.
    vec2 diskStructure(float r, float phi, float omega) {
      if (uDiskDetail <= 0.0) return vec2(1.0);
      float l = log(r / uDiskIn);
      float rings = 0.6 * vnoise(vec3(9.0 * l, 3.1, 0.5)) + 0.4 * vnoise(vec3(23.0 * l, 7.3, 2.5));
      // The layers' deviations from the mean are blended with weights normalized in quadrature, which keeps the
      // contrast of the sum steady through a cross-fade.
      vec2 acc = vec2(0.0);
      float w2 = 0.0;
      for (int k = 0; k < 2; k++) {
        vec4 L = uDiskLayer[k];
        if (L.x <= 0.0) continue;
        float psi = phi - omega * max(uDiskAge + L.y, 0.0);
        vec2 cs = vec2(cos(psi), sin(psi));
        vec3 seed = L.z * vec3(37.1, 11.3, 7.7);
        // Streaks: about 10 per turn along the orbit and 20 per e-fold of radius across it, in four octaves.
        vec3 q = vec3(1.6 * cs, 22.0 * l) + seed;
        float s = 0.0, amp = 0.5;
        for (int o = 0; o < 4; o++) { s += amp * vnoise(q); q = q * 2.03 + 1.7; amp *= 0.5; }
        acc += L.x * vec2(s / 0.9375 - 0.5, vnoise(vec3(3.0 * cs, 6.0 * l) + seed + 5.0) - 0.5);
        w2 += L.x * L.x;
      }
      acc /= sqrt(max(w2, 1e-6));
      float cover = smoothstep(0.18, 0.82, 0.6 * (acc.x + 0.5) + 0.4 * rings);
      return vec2(mix(1.0, 0.12 + 0.88 * cover, uDiskDetail), 1.0 + 0.3 * uDiskDetail * acc.y);
    }
    // A disk crossing: (rgb, coverage). The gas orbits on circular geodesics; the observed specific intensity is
    // that of a blackbody at g T (I_nu / nu^3 is invariant), g the ratio of observed to emitted frequency.
    vec4 disk(vec3 h, vec4 p, float puObs) {
      float r = blR(h);
      float x = uDiskIn / r;
      // Thin-disk temperature profile r^-3/4 (1 - sqrt(r_in/r))^1/4, peak normalized to 1.
      float temp = pow(x, 0.75) * pow(max(1.0 - sqrt(x), 0.0), 0.25) / 0.488;
      float M = 0.5 * uRs;
      float Mr = max(M * r - uQ2, 0.0);
      float omega = sqrt(Mr) / (r * r + uA * sqrt(Mr));
      float g = 1.0;
      if (uShift == 1) {
        vec4 v = vec4(1.0, -omega * h.y, omega * h.x, 0.0);
        float f; vec4 k; ks(h, f, k);
        float kv = dot(k, v);
        float gvv = -1.0 + dot(v.yzw, v.yzw) + f * kv * kv;
        if (gvv < 0.0) {
          vec4 u = v / sqrt(-gvv);
          g = puObs / dot(p, u);
        }
      }
      // Kerr-Schild x + iy = (r + ia) e^(i phi) in the equatorial plane.
      float phi = atan(h.y, h.x) - atan(uA, r);
      vec2 st = diskStructure(r, phi, omega);
      float cover = st.x * (1.0 - smoothstep(0.8 * uDiskOut, uDiskOut, r));
      return vec4(planck(g * uDiskTPeak * temp * st.y) * (0.4 * cover), cover);
    }

    // --- Gas as a volume (uDebris) ---
    // Angular velocity of circular orbits about z at cylindrical radius rc (equatorial Kerr, prograde).
    float orbitOmega(float rc, float M, float a) { return sqrt(M) / (pow(rc, 1.5) + a * sqrt(M)); }
    // Structure of the gas: a two-armed spiral (weight spiralW; a merger's gas has its own) and streaks along
    // the orbits, both sheared by the differential rotation (the disk's two cross-fading layers, so each radius
    // turns at its own rate with the winding bounded), the streaks varying a little with height zn = z / H; and
    // rings. Large-scale structure survives the long sight lines through the gas; streaks mostly average out.
    // (density, temperature factor).
    vec2 gasStructure(float rc, float phi, float zn, float omega, float spiralW) {
      if (uDiskDetail <= 0.0) return vec2(1.0);
      float l = log(rc / uDebIn);
      float rings = vnoise(vec3(9.0 * l, 3.1 + 0.4 * zn, 0.5));
      float acc = 0.0, spiral = 0.0, w2 = 0.0;
      for (int k = 0; k < 2; k++) {
        vec4 L = uDiskLayer[k];
        if (L.x <= 0.0) continue;
        float psi = phi - omega * max(uDiskAge + L.y, 0.0);
        vec3 q = vec3(1.6 * cos(psi) + 0.3 * zn, 1.6 * sin(psi) - 0.2 * zn, 14.0 * l) + L.z * vec3(37.1, 11.3, 7.7);
        acc += L.x * (0.65 * vnoise(q) + 0.35 * vnoise(2.03 * q + 1.7) - 0.5);
        spiral += L.x * (0.5 * cos(2.0 * psi + 2.0 * l + 2.3 * L.z) + 0.25 * cos(psi - 1.7 * L.z));
        w2 += L.x * L.x;
      }
      acc /= sqrt(max(w2, 1e-6));
      spiral *= spiralW / sqrt(max(w2, 1e-6));
      float dens = smoothstep(0.1, 0.9, 0.5 + 0.55 * spiral + 0.8 * acc + 0.3 * (rings - 0.5));
      return vec2(mix(1.0, 0.05 + 1.6 * dens, uDiskDetail), 1.0 + uDiskDetail * (0.35 * spiral + 0.6 * acc + 0.2 * (rings - 0.5)));
    }
    // Optical depth per unit length of the gas at x, its temperature factor, and tint: the colour of the stars
    // the young gas came from (rgb) and the fraction of the gas it makes up (a).
    // A merger's gas: the particles' density from uGasTex, with finer streaks on top (half strength) that turn
    // with the disk's layers. A thick disk: a Gaussian in height with a vertical optical depth of 3 r_in / r,
    // opaque only near the hole, so elsewhere its brightness follows the density of its structure.
    // Either way the midplane is up to 1.25 times hotter than the cooler atmosphere.
    float gas(vec3 x, out float tf, out vec4 tint) {
      tf = 1.0;
      tint = vec4(0.0);
      float rc = length(x.xy);
      float phi = atan(x.y, x.x);
      float H = uDebH * max(rc, 1e-3);
      float vert = exp(-0.5 * x.z * x.z / (H * H));
      if (uDebris == 1) {
        vec3 tc = (x - uGasMin) / uGasSize;
        if (any(lessThan(tc, vec3(0.0))) || any(greaterThan(tc, vec3(1.0)))) return 0.0;
        vec4 g = texture(uGasTex, tc);
        if (g.r <= 1e-6) return 0.0;
        vec3 edge = min(tc, 1.0 - tc);
        // (Faded beyond the disk's outer radius, where only a few particles have wandered.)
        float fade = smoothstep(0.0, 0.06, min(min(edge.x, edge.y), edge.z)) * (1.0 - smoothstep(0.8 * uDebOut, 1.2 * uDebOut, length(x)));
        vec2 st = gasStructure(max(rc, 0.5 * uDebIn), phi, x.z / H, orbitOmega(max(rc, 0.5 * uDebIn), uDebM, uDebA), 0.0);
        tf = mix(1.0, st.y, 0.5) * (0.8 + 0.45 * vert);
        float young = g.g + g.b;
        if (young > 0.0) tint = vec4((g.g * uStarColour1 + g.b * uStarColour2) / young, min(young / g.r, 1.0));
        return uGasScale * g.r * fade * mix(1.0, st.x, 0.5);
      }
      if (uFill <= 0.0 || rc < 0.8 * uDebIn || rc > uDebOut || abs(x.z) > 3.5 * H) return 0.0;
      float radial = smoothstep(0.8 * uDebIn, 1.5 * uDebIn, rc) * (1.0 - smoothstep(0.8 * uFront, uFront, rc))
        * (1.0 - smoothstep(0.6 * uDebOut, uDebOut, rc)) * (uDebIn / rc);
      if (radial <= 1e-3) return 0.0;
      vec2 st = gasStructure(rc, phi, x.z / H, orbitOmega(rc, uDebM, uDebA), 1.0);
      tf = st.y * (0.8 + 0.45 * vert);   // a hotter midplane under a cooler atmosphere
      return uFill * radial * st.x * vert * (3.0 / (2.5066 * H));
    }
    // The gas's source function: a blackbody at the thin-disk temperature of its radius (inner edge uDebIn),
    // seen at g T for gas on circular orbits, scaled a little below an opaque thin disk's (a thick torus shows
    // much more glowing surface).
    vec3 gasSource(vec3 x, vec4 p, float puObs, float tf) {
      float rc = max(length(x.xy), uDebIn);
      float xr = uDebIn / rc;
      float temp = pow(xr, 0.75) * pow(max(1.0 - sqrt(xr), 0.0), 0.25) / 0.488;
      float g = 1.0;
      if (uShift == 1) {
        float om = orbitOmega(rc, uDebM, uDebA);
        vec4 v = vec4(1.0, -om * x.y, om * x.x, 0.0);
        float gvv = dot(v, metricAt(x) * v);
        if (gvv < 0.0) g = puObs / dot(p, v / sqrt(-gvv));
      }
      return planck(g * uDiskTPeak * temp * tf) * 0.22;
    }

    // Whether the straight segment from P to P + D passes the hole (closest approach inside uShadowR, softened).
    float clearOfHole(vec3 P, vec3 D) {
      if (uShadowR <= 0.0) return 1.0;
      float t = clamp(-dot(P, D) / max(dot(D, D), 1e-12), 0.0, 1.0);
      return smoothstep(0.9 * uShadowR, 1.1 * uShadowR, length(P + t * D));
    }
    // Light on a surface point P with normal n from the accretion disk, treated as 6 rings of 16 glowing points
    // (each one an area of the disk, with the disk's blackbody radiance at its radius and the mean coverage of its
    // structure) along straight paths, blocked by the hole. The points emit in all directions (a thin disk seen
    // edge-on would leave a body in its plane dark), and the distance is softened by a third of the radius, so a
    // body inside the disk is lit smoothly.
    vec3 diskIrradiance(vec3 P, vec3 n) {
      vec3 E = vec3(0.0);
      float l0 = log(uDiskIn), l1 = log(uDiskOut);
      for (int i = 0; i < 6; i++) {
        float r = exp(mix(l0, l1, (float(i) + 0.5) / 6.0));
        float x = uDiskIn / r;
        float temp = pow(x, 0.75) * pow(max(1.0 - sqrt(x), 0.0), 0.25) / 0.488;
        float dA = r * r * (l1 - l0) / 6.0 * (2.0 * PI / 16.0);
        vec3 L = planck(uDiskTPeak * temp) * (0.4 * 0.6 * (1.0 - smoothstep(0.8 * uDiskOut, uDiskOut, r)) * dA);
        for (int j = 0; j < 16; j++) {
          float ph = (float(j) + 0.5 * float(i % 2)) * (2.0 * PI / 16.0);
          vec3 D = vec3(r * cos(ph) - uA * sin(ph), r * sin(ph) + uA * cos(ph), 0.0) - P;   // (r + ia) e^(i phi)
          float d2 = dot(D, D);
          float c = dot(n, D) * inversesqrt(d2);
          if (c <= 0.0) continue;
          E += L * (c * clearOfHole(P, D) / (d2 + 0.11 * r * r));
        }
      }
      return E;
    }
    // A distant star's light reaching P: eclipsed by the hole, dimmed by the disk's coverage where the path crosses it.
    float starLight(vec3 P) {
      float vis = clearOfHole(P, uSunDir * 1e6);
      if (uDisk == 1 && P.z * uSunDir.z < 0.0) {
        vec3 Q = P - uSunDir * (P.z / uSunDir.z);
        float r = length(Q.xy);
        if (r >= uDiskIn && r <= uDiskOut) {
          float Mr = max(0.5 * uRs * r - uQ2, 0.0);
          float omega = sqrt(Mr) / (r * r + uA * sqrt(Mr));
          float cover = diskStructure(r, atan(Q.y, Q.x) - atan(uA, r), omega).x * (1.0 - smoothstep(0.8 * uDiskOut, uDiskOut, r));
          vis *= 1.0 - 0.88 * cover;
        }
      }
      return vis;
    }
    // What a body's surface at P (normal n; toCam towards the camera) receives; its colour is the albedo times
    // this. Not ray traced: straight light paths. The disk's light is tone-mapped (1 - exp(-E)) so a body stays
    // visible far from it while keeping the light's colour and direction.
    vec3 bodyLight(vec3 P, vec3 n, vec3 toCam) {
      if (uLightMode == 2) return vec3((0.35 + 0.65 * max(dot(n, toCam), 0.0)) * 1.4);
      vec3 light;
      if (uLightMode == 0 && uDisk == 1) {
        vec3 E = diskIrradiance(P, n);
        float I = dot(E, vec3(0.2126, 0.7152, 0.0722));
        light = I > 1e-9 ? E * (1.5 * (1.0 - exp(-I)) / I) : vec3(0.0);
      } else {
        float c = max(dot(n, uSunDir), 0.0);
        light = vec3(1.0, 0.96, 0.9) * (1.4 * c * (c > 0.0 ? starLight(P) : 0.0));
      }
      return light + uBodyGlow;
    }
    // The Roche potential and its gradient (the outward normal of an equipotential surface).
    float roche(vec3 x) {
      vec3 c = x - uCM;
      return -uM1 / length(x - uX1) - uM2 / length(x - uX2) - 0.5 * uOmega2 * dot(c.xy, c.xy);
    }
    vec3 rocheGrad(vec3 x) {
      vec3 d1 = x - uX1, d2 = x - uX2, c = x - uCM;
      float l1 = length(d1), l2 = length(d2);
      return uM1 * d1 / (l1 * l1 * l1) + uM2 * d2 / (l2 * l2 * l2) - uOmega2 * vec3(c.xy, 0.0);
    }
    // Spherical stars (uTides = 0): which one x is inside, -1 none.
    int insideSphere(vec3 x) {
      return length(x - uX1) < uR1 ? 0 : (length(x - uX2) < uR2 ? 1 : -1);
    }
    // The tidal field D (uTides = 1). Inside a bounding sphere, the first-order distance (phi - level) / |grad phi|
    // to the equipotential surface (the gradient floored, as it vanishes at the Lagrange points); outside both,
    // the distance to them. Every parameter changes continuously in time, so the surface does too, through the
    // touching of the stars (the levels reach the inner Lagrange point's) and the blend into the ball.
    float shapeField(vec3 x) {
      float dBall = blR(x) - uBall;
      if (uBlend >= 0.999) return dBall;
      float dR = min(length(x - uX1) - uB1, length(x - uX2) - uB2);
      if (dR <= 0.0) {
        float level = dot(x - uL1, uX1 - uX2) > 0.0 ? uLevel1 : uLevel2;
        dR = (roche(x) - level) / max(length(rocheGrad(x)), uGFloor);
      }
      return uBlend > 0.0 ? mix(dR, dBall, uBlend) : dR;
    }
    vec3 shapeNormal(vec3 x) {
      float h = 2e-3 * max(max(uB1, uB2), uBall) + 1e-5;
      const vec2 e = vec2(1.0, -1.0);
      vec3 n = e.xyy * shapeField(x + h * e.xyy) + e.yyx * shapeField(x + h * e.yyx) + e.yxy * shapeField(x + h * e.yxy) + e.xxx * shapeField(x + h * e.xxx);
      return length(n) > 0.0 ? normalize(n) : normalize(x);
    }
    // Colour of the surface at x: each star's on its side of the inner Lagrange plane (blended across the neck),
    // whitened by the heat, and dimmed by the gravitational redshift (sqrt(-g_tt) at x relative to the camera's
    // gttCam, cubed, relative to 0.7, about a neutron star's own): the remnant fades to black as it collapses
    // onto the horizon, where static emitters cease to exist.
    vec3 shapeColour(vec3 x, float gttCam) {
      vec3 axis = uX1 - uX2;
      float s = length(axis) > 1e-3 ? clamp(0.5 + dot(x - uL1, normalize(axis)) / 1.5, 0.0, 1.0) : 0.5;
      vec3 c = mix(mix(uStarColour2, uStarColour1, s), vec3(1.6, 1.5, 1.35), uHeat);
      float gtt = metricAt(x)[0][0];
      float red = gtt < 0.0 ? clamp(sqrt(gtt / gttCam) / 0.7, 0.0, 1.0) : 0.0;
      return c * red * red * red;
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
      bool done = false, hit = false, escaped = false;
      bool tube = uMode == 3 && inTube(x);
      vec3 v = velocity(x, p);
      vec3 xPrev = x;
      for (int i = 0; i < 4096; i++) {
        if (i >= uMaxSteps) break;
        float r = blR(x);
        // Step size from the distance to the nearest horizon or surface, with a floor proportional to the
        // distance from that body's centre (rays winding onto a horizon still make progress).
        float rNear = r - uHorizon, rCentre = r, dsStar = 0.0;
        if (uBinary == 1 || uRemnant == 1) {
          // Each body: a horizon at r_a = 2 m_a (as for one Kerr-Schild hole), or a star: a sphere, or (uTides) the
          // tidal field, which after the merger is the remnant, also once the Kerr metric has taken over (uRemnant).
          float r1 = length(x - uX1), r2 = length(x - uX2);
          bool stars = uR1 > 0.0 || uR2 > 0.0;
          bool tidal = uTides == 1;
          float h1 = 2.0 * uM1 * 1.01, h2 = 2.0 * uM2 * 1.01;
          float D = tidal ? shapeField(x) : 1e30;
          bool star = tidal ? D < 0.0 : stars && insideSphere(x) >= 0;
          bool hole = uBinary == 1 && ((uR1 == 0.0 && r1 < h1) || (uR2 == 0.0 && r2 < h2));
          if (star || hole) {
            if (star) {
              // The surface by bisection along the last step, its normal and colour.
              vec3 a = xPrev, b = x;
              for (int k = 0; k < 6; k++) { vec3 m = 0.5 * (a + b); if (tidal ? shapeField(m) < 0.0 : insideSphere(m) >= 0) b = m; else a = m; }
              vec3 n, c;
              if (tidal) { n = shapeNormal(b); c = shapeColour(b, A); }
              else { int s = insideSphere(b); n = normalize(b - (s == 0 ? uX1 : uX2)); c = s == 0 ? uStarColour1 : uStarColour2; }
              colour += transmit * c * (0.35 + 0.65 * max(dot(n, -normalize(v)), 0.0));
            }
            done = true; hit = true; break;
          }
          // Steps. Horizons need steps that shrink with the distance to them (rays wind onto them); star surfaces
          // do not, and rays grazing them would run out of steps: half the distance to the surface, with a floor,
          // under the curvature control below.
          float dHole = min(uR1 == 0.0 ? r1 - h1 : 1e30, uR2 == 0.0 ? r2 - h2 : 1e30);
          float dS = tidal ? D : min(uR1 > 0.0 ? r1 - uR1 : 1e30, uR2 > 0.0 ? r2 - uR2 : 1e30);
          if (stars) dsStar = max(0.5 * dS, 0.01 * max(max(uB1, uB2), max(max(uR1, uR2), uBall)));
          if (uBinary == 1) { rNear = dHole; rCentre = min(r1, r2); }
        } else if (uRrem > 0.0 && r < uRrem) {
          vec3 n = normalize(x);
          colour += transmit * uStarColourRem * (0.35 + 0.65 * max(dot(n, -normalize(v)), 0.0));
          done = true; hit = true; break;
        } else if (uRrem > 0.0) rNear = r - uRrem;
        if (uHorizon > 0.0 && r < uHorizon * 1.01) { done = true; break; }
        if (r > uEscape && dot(x, v) > 0.0) { escaped = true; done = true; break; }
        // Outside every body's reach and the photon orbits, an outgoing ray never comes back to a body.
        if (uMode >= 2 && r > uBodyReach && dot(x, v) > 0.0) { done = true; break; }
        // The step: a fraction of the distance to the nearest horizon, or (curvature control) of the distance to
        // the nearest mass, with a floor; star surfaces (dsStar) can only shorten it.
        float ds = uStepFrac * max(min(rNear, rCentre), 0.02 * rCentre + 1e-4);
        if (dsStar > 0.0) ds = min(ds, dsStar);
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
        // Gas as a volume: emission and absorption at one point of the step, placed at random along it (the
        // refinement samples average the placements).
        if (uDebris > 0) {
          float tm = 0.05 + 0.9 * hash13(vec3(gl_FragCoord.xy + 1024.0 * uJitter, float(i) + uFrameSeed));
          vec3 xm = mix(x, xn, tm);
          float tf; vec4 tint;
          float rho = gas(xm, tf, tint);
          if (rho > 1e-4) {
            float absorbed = 1.0 - exp(-rho * segLen);
            if (uMode != 3) colour += transmit * absorbed * mix(gasSource(xm, mix(p, pn, tm), puObs, tf), 0.8 * tint.rgb, tint.a);
            transmit *= 1.0 - absorbed;
          }
        }
        // Thin disk in the equatorial plane z = 0.
        float diskT = 1e30;
        if (uDisk == 1 && uDebris == 0 && x.z * xn.z < 0.0) {
          float t = x.z / (x.z - xn.z);
          vec3 hp = mix(x, xn, t);
          float rh = blR(hp);
          if (rh >= uDiskIn && rh <= uDiskOut) diskT = t * segLen;
        }
        if (diskT < hitT) {
          float t = diskT / max(segLen, 1e-12);
          vec4 dc = disk(mix(x, xn, t), mix(p, pn, t), puObs);
          if (uMode != 3) colour += transmit * dc.rgb;
          transmit *= 1.0 - 0.88 * dc.a;
        }
        if (hitBody >= 0) {
          vec3 hp = x + sd * hitT;
          vec3 n = normalize(hp - uBody[hitBody].xyz);
          colour += transmit * uBodyColor[hitBody] * bodyLight(hp, n, -sd);
          done = true; hit = true;
          break;
        }
        xPrev = x; x = xn; p = pn;
        v = velocity(x, p);
        if (transmit < 0.02) { done = true; break; }
      }
      // Out of steps: a ray still near the hole is winding onto the horizon (traced backward in time it
      // approaches the past horizon only asymptotically in ingoing Kerr-Schild coordinates), so it is
      // captured; a ray still far out is given the sky in its current direction.
      bool nearHole = uBinary == 1
        ? (uR1 == 0.0 && length(x - uX1) < 8.0 * uM1) || (uR2 == 0.0 && length(x - uX2) < 8.0 * uM2)
        : blR(x) <= 4.0 * max(uHorizon, 1e-6);
      // The sky: the ray's final direction, and its change from pixel to pixel (screen-space derivatives, taken here
      // in uniform control flow), which is each pixel's footprint on the sky.
      bool skyShown = escaped || (!done && !nearHole);
      vec3 skyDir = normalize(v);
      vec3 sx = dFdx(skyDir), sy = dFdy(skyDir);
      float camArea = length(cross(dFdx(d), dFdy(d)));
      float shown = skyShown ? 1.0 : 0.0;
      bool footprintOk = abs(dFdx(shown)) + abs(dFdy(shown)) == 0.0;
      if (skyShown) colour += transmit * sky(skyDir, sx, sy, camArea, footprintOk);
      // Linear radiance with alpha 1: samples are summed in a float buffer and averaged by alpha at display.
      if (uMode == 2) gl_FragColor = hit ? vec4(colour, 1.0) : vec4(0.0);
      else gl_FragColor = vec4(colour, tube ? 2.0 : 1.0);
    }`;

  window.GRStyles = window.GRStyles || [];
  window.GRStyles.push({
    id: 'raytrace',
    label: 'Ray tracer',
    order: 5,
    description: 'Null geodesics traced per pixel on the GPU: the shadow and photon ring, the lensed sky, an accretion disk (thin, or a volume of gas) glowing as a blackbody (its temperature set by the hole\'s mass and accretion rate) with streaks and rings and with Doppler and gravitational shifts, the bodies as a camera sees them, and for merging neutron stars their tidal shapes, shed matter and debris torus. Drag to orbit; or ride along with a body.',
    options: [
      { id: 'resolution', label: 'max resolution', type: 'range', min: 0.25, max: 1, step: 0.05, value: 1, help: 'Fraction of the display\'s device pixels traced once the view is still' },
      { id: 'samples', label: 'refine samples', type: 'range', min: 1, max: 64, step: 1, value: 24, help: 'Jittered full-resolution samples averaged while the view is still (anti-aliasing)' },
      { id: 'steps', label: 'max steps', type: 'range', min: 100, max: 2000, step: 50, value: 600 },
      { id: 'accuracy', label: 'step size', type: 'range', min: 0.01, max: 0.12, step: 0.005, value: 0.04, help: 'Spatial step as a fraction of the distance to the horizon' },
      // No disk by default without a mass to orbit (Minkowski scenes such as the twin paradox).
      // A disk around a mass, or the debris of a neutron-star merger that leaves a black hole.
      { id: 'disk', label: 'accretion disk', type: 'checkbox', value: (data) => { const f = familyOf(data); return f && (f.rs > 0 || (f.binary && f.binary.R1 > 0 && f.binary.R2 > 0 && f.binary.Rrem === 0)) ? 1 : 0; } },
      { id: 'diskOuter', label: 'disk outer (rs)', type: 'range', min: 3, max: 30, step: 0.5, value: (data) => (neutronStars(data) ? 14 : 8), show: diskShown,
        help: 'Outer radius of the disk, in units of the Schwarzschild radius (of M for a binary)' },
      { id: 'diskThick', label: 'disk thickness (H/r)', type: 'range', min: 0, max: 0.5, step: 0.02, value: (data) => (neutronStars(data) ? 0.22 : 0), show: diskShown,
        help: 'Scale height over radius: 0 is a thin disk; above it the disk is a volume of glowing, absorbing gas (slower to trace). The debris torus of a neutron-star merger is always a volume' },
      { id: 'diskColour', label: 'disk temperature', type: 'select', value: 0, show: diskShown,
        choices: [{ value: 0, label: 'from mass and accretion rate' }, { value: 1, label: 'set the peak' }],
        help: 'Each ring of the disk glows as a blackbody at the thin-disk temperature for its radius; the scale comes from the hole\'s mass and accretion rate, or is set directly' },
      // A scene with a physical mass (Examples.BinaryMerger's gwSource) starts with it; otherwise M87*'s 6.3e9.
      { id: 'diskMass', label: 'hole mass (log10 M☉)', type: 'range', min: 0, max: 10.5, step: 0.1, show: (d, o) => diskShown(d, o) && !(o && o.diskColour == 1),
        value: (data) => { const g = data && data.records && (data.records.gwSource || [])[0]; return g && g.totalMassSolar > 0 ? Math.round(10 * Math.log10(g.totalMassSolar)) / 10 : 9.8; } },
      { id: 'diskRate', label: 'accretion (log10 Eddington)', type: 'range', min: -6, max: 0, step: 0.1, value: -4.5, show: (d, o) => diskShown(d, o) && !(o && o.diskColour == 1),
        help: 'Mass accretion rate as a fraction of the Eddington rate (log10)' },
      { id: 'diskTemp', label: 'peak temperature (K)', type: 'range', min: 1000, max: 40000, step: 100, value: 5500, show: (d, o) => diskShown(d, o) && !!(o && o.diskColour == 1) },
      { id: 'diskDetail', label: 'disk structure', type: 'range', min: 0, max: 1, step: 0.05, value: 0.8, show: diskShown,
        help: 'Streaks along the orbits, rings and gaps, warm and cool patches; 0 is a smooth disk' },
      { id: 'diskMotion', label: 'disk rotates', type: 'checkbox', value: 0, show: diskShown,
        help: 'The pattern turns at the orbital speed of each radius with the clock. Every clock change then re-traces the whole image, so this is meant for rendered videos (the Render dialog has the same switch)' },
      { id: 'tides', label: 'tides and debris', type: 'checkbox', value: 1, show: (d) => neutronStars(d),
        help: 'Neutron stars: their tidal (Roche) shapes as they close in, the arms of matter they shed, the merged remnant collapsing into the hole, and the debris torus' },
      { id: 'startSep', label: 'start at separation (M)', type: 'range', min: 4, max: 30, step: 0.5, value: (d) => (neutronStars(d) ? 12 : 30), show: (d) => mergerScene(d),
        help: 'The animation (and the Render dialog) starts when the bodies are this far apart; at or above the first recorded separation it starts at the beginning' },
      { id: 'shift', label: 'Doppler + redshift', type: 'checkbox', value: 1 },
      { id: 'sky', label: 'sky', type: 'select', choices: [{ value: 0, label: 'stars' }, { value: 1, label: 'stars + grid' }, { value: 2, label: 'checker' }], value: 1 },
      { id: 'bodySize', label: 'body size (rs)', type: 'range', min: 0.05, max: 1.5, step: 0.05, value: 0.25 },
      { id: 'bodyLight', label: 'bodies lit by', type: 'select', value: 0,
        choices: [{ value: 0, label: 'the accretion disk' }, { value: 1, label: 'a distant star' }, { value: 2, label: 'the camera (flat)' }],
        help: 'Light on the bodies along straight paths (not lensed): the disk\'s own blackbody glow, blocked by the hole, or a distant star that the hole can eclipse and the disk can shade. Without a disk the star lights them' },
      { id: 'sunAz', label: 'star azimuth (deg)', type: 'range', min: 0, max: 360, step: 5, value: 100, show: (d, o) => !o || o.bodyLight == 1 || (o.bodyLight == 0 && !o.disk) },
      { id: 'sunEl', label: 'star elevation (deg)', type: 'range', min: -90, max: 90, step: 5, value: 30, show: (d, o) => !o || o.bodyLight == 1 || (o.bodyLight == 0 && !o.disk) },
      { id: 'bodyGlow', label: 'night-side glow', type: 'range', min: 0, max: 0.5, step: 0.01, value: 0.06, help: 'A faint glow of the bodies themselves (not physical), so their dark side stays visible' },
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
      // The Roche potential of point masses m1, m2 at X1, X2 in the frame rotating at sqrt(omega2) about their
      // centre of mass (orbital plane z = 0). Each star's own surface level is the potential at its radius R on the
      // side away from the companion; q is how far the higher of them is from the inner Lagrange point's, relative
      // to it (q = 0: it reaches L1). From q = -0.08 to -0.01 the two levels blend into their mean (opts.c forces
      // the weight), so the stars reach L1 at one level and the neck opens continuously; the levels stay below the
      // outer Lagrange points', so the surfaces are closed. Each bounding sphere holds its star's part of the
      // surface (its reach along 13 directions and to the rim of the neck, with a margin), short of the outer
      // Lagrange point on its side; rEq is the radius of a sphere of the pair's size.
      function roche(X1, X2, m1, m2, R1, R2, omega2, opts = {}) {
        const cm = X1.clone().multiplyScalar(m1).add(X2.clone().multiplyScalar(m2)).divideScalar(m1 + m2);
        const a = Math.max(X1.distanceTo(X2), 1e-6), e = X2.clone().sub(X1).divideScalar(a);
        const P = new THREE.Vector3();
        const phiAt = (V) => { const cx = V.x - cm.x, cy = V.y - cm.y; return -m1 / Math.max(V.distanceTo(X1), 1e-9) - m2 / Math.max(V.distanceTo(X2), 1e-9) - 0.5 * omega2 * (cx * cx + cy * cy); };
        const phi = (s) => phiAt(P.copy(X1).addScaledVector(e, s));
        const argmax = (lo, hi) => { for (let k = 0; k < 50; k++) { const u = lo + (hi - lo) / 3, w = hi - (hi - lo) / 3; if (phi(u) < phi(w)) lo = u; else hi = w; } return 0.5 * (lo + hi); };
        const sL1 = argmax(1e-3 * a, a * (1 - 1e-3)), sL3 = argmax(-4 * a - 4 * R1, -1e-3 * a), sL2 = argmax(a * (1 + 1e-3), 5 * a + 4 * R2);
        const phiL1 = phi(sL1), phiOuter = Math.min(phi(sL2), phi(sL3)), cap = phiOuter - 0.02 * Math.abs(phiOuter);
        const own1 = phi(-R1), own2 = phi(a + R2);
        const q = (Math.max(own1, own2) - phiL1) / Math.abs(phiL1);
        const c = opts.c ?? smooth(q, -0.08, -0.01);
        const mean = 0.5 * (own1 + own2);
        const level1 = Math.min(own1 + c * (mean - own1), cap), level2 = Math.min(own2 + c * (mean - own2), cap);
        const L1 = X1.clone().addScaledVector(e, sL1);
        const res = { cm, L1, level1, level2, q, connected: Math.max(level1, level2) > phiL1 };
        if (opts.quick) return res;
        const perp = new THREE.Vector3(-e.y, e.x, 0).normalize(), up = new THREE.Vector3(0, 0, 1);
        // Reach of star i's part of the surface from its centre X along u: the first point at or above its level,
        // or the inner Lagrange plane (at distance toPlane along f, the direction to the companion).
        const V = new THREE.Vector3();
        const reach = (X, u, level, f, toPlane, tMax) => {
          const uf = u.dot(f), tPlane = uf > 1e-6 ? toPlane / uf : Infinity, end = Math.min(tPlane, tMax);
          const n = 60, h = end / n;
          for (let k = 1; k <= n; k++) {
            if (phiAt(V.copy(X).addScaledVector(u, k * h)) >= level) {
              let lo = (k - 1) * h, hi = k * h;
              for (let j = 0; j < 30; j++) { const m = 0.5 * (lo + hi); if (phiAt(V.copy(X).addScaledVector(u, m)) >= level) hi = m; else lo = m; }
              return hi;
            }
          }
          return end;
        };
        const ext = [];
        for (const [X, f, level, toPlane, outer] of [[X1, e, level1, sL1, -sL3], [X2, e.clone().negate(), level2, a - sL1, sL2 - a]]) {
          const b = f.clone().negate();
          const dirs = [b, f, perp, perp.clone().negate(), up];
          for (const [p1, p2] of [[b, perp], [b, perp.clone().negate()], [b, up], [f, perp], [f, perp.clone().negate()], [f, up], [perp, up], [perp.clone().negate(), up]]) dirs.push(p1.clone().add(p2).normalize());
          let far = 0;
          const r = {};
          dirs.forEach((u, k) => { const d = reach(X, u, level, f, toPlane, outer); far = Math.max(far, d); if (k === 0) r.back = d; if (k === 2 || k === 3) r.perp = Math.max(r.perp || 0, d); if (k === 4) r.z = d; });
          if (level > phiL1) {
            // The neck: its rim in the Lagrange plane, along the in-plane perpendicular and z.
            for (const u of [perp, up]) {
              let lo = 0, hi = 3 * a + 3 * Math.max(R1, R2);
              for (let j = 0; j < 40; j++) { const m = 0.5 * (lo + hi); if (phiAt(V.copy(L1).addScaledVector(u, m)) >= level) hi = m; else lo = m; }
              far = Math.max(far, Math.hypot(toPlane, hi));
            }
          }
          r.B = Math.min(1.15 * far + 0.05, 0.97 * outer);
          ext.push(r);
        }
        res.B1 = ext[0].B; res.B2 = ext[1].B;
        res.rEq = Math.cbrt(0.5 * (ext[0].back + a + ext[1].back) * Math.max(ext[0].perp, ext[1].perp) * Math.max(ext[0].z, ext[1].z));
        return res;
      }
      // A merger in time (ct, in M), from the recorded orbit: the separation and orientation of the pair before
      // tMerge, and their values at it (extrapolated from the last samples). For neutron stars (option "tides and
      // debris"), the sketch's timeline:
      //   - the tidal shapes (roche above), with the orbit's own angular velocity for the rotating frame;
      //   - the outer tips shed matter from tShed, when the surfaces come within 12% of the inner Lagrange point's
      //     potential; the arms are launched until tLaunchEnd;
      //   - after the merger the pair, still turning at the orbit's rate, closes up over tJoin (0.6 orbits) as a
      //     smooth (Hermite) continuation of the plunge, its shape blending into a ball, while the remnant's
      //     Kerr term fades into the two-centre metric (binaryMetric), so the lensed sky slows from turning with
      //     the pair to the spinning hole's fixed pattern without a stop in between; at tJoin the metric is the
      //     remnant's Kerr metric;
      //   - the ball shrinks onto the horizon, dimmed by its redshift, by tEnd (0.9 orbits later);
      //   - the torus fills in over the next orbits while the arms wind up; shed gas keeps its star's colour at
      //     first and turns to the gas's blackbody colour over 0.8 orbits.
      // slow: the span the Render dialog can slow down, from shedding to the end of the collapse (for black holes
      // the last 1.5 orbits and 100 M of ringdown).
      function buildStory() {
        const tm = bin.tMerge, Mt = bin.m1 + bin.m2;
        if (!(tm < 1e29) || data.trackSpecs.length < 2) return null;
        const s1 = data.trackSpecs[0].samples, s2 = data.trackSpecs[1].samples;
        const T = [], A = [], PH = [];
        for (let k = 0, n = Math.min(s1.length, s2.length); k < n; k++) {
          const t = s1[k].x[0];
          if (t >= tm) break;
          const d = samplePhys(s1[k].x).sub(samplePhys(s2[k].x));
          let ph = Math.atan2(d.y, d.x);
          if (PH.length) { const p0 = PH[PH.length - 1]; ph = p0 + Math.atan2(Math.sin(ph - p0), Math.cos(ph - p0)); }
          T.push(t); A.push(d.length()); PH.push(ph);
        }
        const L = T.length - 1;
        if (L < 4) return null;
        const OM = PH.map((_, k) => (PH[Math.min(k + 1, L)] - PH[Math.max(k - 1, 0)]) / (T[Math.min(k + 1, L)] - T[Math.max(k - 1, 0)]));
        const at = (arr, t) => {
          if (t >= T[L]) return arr[L] + (arr[L] - arr[L - 1]) / (T[L] - T[L - 1]) * (t - T[L]);
          if (t <= T[0]) return arr[0];
          let lo = 0, hi = L;
          while (hi - lo > 1) { const m = (lo + hi) >> 1; if (T[m] <= t) lo = m; else hi = m; }
          return arr[lo] + (arr[hi] - arr[lo]) * (t - T[lo]) / (T[hi] - T[lo]);
        };
        const am = Math.max(at(A, tm), 0.5 * A[L]), phm = at(PH, tm), omega = at(OM, tm);
        const adot = (A[L] - A[L - 1]) / (T[L] - T[L - 1]);
        const P = 2 * Math.PI / Math.max(Math.abs(omega), 1e-9);
        const st = { tm, T, A, at, am, phm, omega, adot, P, omega2At: (t) => at(OM, Math.min(t, tm)) ** 2 };
        st.startAt = (sep) => { if (sep >= A[0]) return T[0]; for (let k = 0; k <= L; k++) if (A[k] <= sep) return T[k]; return T[L]; };
        st.slow = [tm - 1.5 * P, tm + 100 * bin.Mf];
        st.slowFactor = 1;
        if (!(bin.R1 > 0 && bin.R2 > 0)) return st;
        const pair = (sep, omega2, opts) => roche(new THREE.Vector3(sep * bin.m2 / Mt, 0, 0), new THREE.Vector3(-sep * bin.m1 / Mt, 0, 0), bin.m1, bin.m2, bin.R1, bin.R2, omega2, opts);
        let tShed = tm - 0.5 * P;
        for (let k = L; k >= 0 && T[k] > tm - 8 * P; k--) {
          if (pair(A[k], OM[k] ** 2, { quick: true }).q >= -0.12) tShed = T[k]; else break;
        }
        st.qm = pair(am, omega * omega, { quick: true }).q;
        st.heatM = 0.35 * smooth(st.qm, -0.05, 0.15);
        st.tShed = tShed;
        st.tJoin = 0.6 * P;
        st.tEnd = st.tJoin + 0.9 * P;
        st.tLaunchEnd = tm + 0.8 * st.tJoin;
        st.dMin = 0.02;
        st.rJoin = pair(st.dMin, omega * omega, { c: 1 }).rEq;
        // The shedding rate (relative), rising from tShed and ending as the remnant closes up.
        st.rate = (te) => smooth(te, tShed, tShed + 0.6 * (tm - tShed)) * (1 - smooth(te, tm + 0.3 * st.tJoin, st.tLaunchEnd));
        // The two masses at time t until they join (x1, y1, x2, y2, separation; the centre of mass at the
        // origin), else null: the recorded orbit, then the pair closing up (a cubic Hermite continuation of the
        // separation from its value and rate at the merger, floored), still turning at the orbit's rate.
        st.pairAt = (t) => {
          let d, ph;
          if (t < tm) { d = at(A, t); ph = at(PH, t); }
          else if (t < tm + st.tJoin) {
            const u = (t - tm) / st.tJoin, k = st.tJoin * adot / am;
            d = Math.max(am * ((2 * u - 3) * u * u + 1 + k * (u - 1) * (u - 1) * u), st.dMin);
            ph = phm + omega * (t - tm);
          } else return null;
          const c = Math.cos(ph), sn = Math.sin(ph), f1 = d * bin.m2 / Mt, f2 = -d * bin.m1 / Mt;
          return [f1 * c, f1 * sn, f2 * c, f2 * sn, d];
        };
        st.rateAt = (t) => (t < tm ? at(OM, t) : omega);
        // Slow motion through the collapse and the torus's filling (3 orbits), where the gas turns fast.
        st.slow = [tShed - 0.25 * P, tm + Math.max(st.tEnd + 0.3 * P, 3 * P)];
        st.slowFactor = 4;
        return st;
      }
      let story;   // built on first use (needs samplePhys)

      // --- The merger's gas as particles -----------------------------------------------------------------
      // GAS_N test particles leave the stars' outer tips (launch times by story.rate, from tShed to tLaunchEnd;
      // 60% from the lighter star): at the tip, spread over its face, moving with the turning pair plus an
      // outward kick of 0.15 to 0.35 c. They move under the pair's gravity (softened point masses on the recorded
      // orbit, then closing up as story.pairAt), then the remnant's at the centre. After the merger the gas
      // dissipates as shocks would: each particle's radial velocity is damped within about half an orbit while its
      // angular momentum is kept, so the gas settles where that angular momentum puts it, and its vertical
      // motion relaxes to a dispersion that holds the torus at the chosen thickness (an Ornstein-Uhlenbeck
      // process; a little in-plane noise spreads it). Gas that falls into a star, the remnant or inside the
      // hole's innermost stable orbit is removed. So the gas leaving the stars is the gas that forms the disk.
      // Leapfrog steps of GAS_DT, deterministic (seeded per particle and step), with a checkpoint every
      // GAS_EVERY steps so the clock can go back; each frame the particles are spread into the GAS_GRID texture
      // (a separable Gaussian over 5 x 5 x 3 voxels, about 0.5 M wide), drifted to the exact time.
      const GAS_N = 30000, GAS_DT = 0.3, GAS_EVERY = 64, GAS_GRID = [128, 128, 32], GAS_KAPPA = 0.1;
      const hashU = (a, b) => {
        let h = Math.imul(a ^ 0x9e3779b9, 0x85ebca6b) ^ Math.imul((b + 0x632be5ab) | 0, 0xc2b2ae35);
        h ^= h >>> 16; h = Math.imul(h, 0x7feb352d); h ^= h >>> 15; h = Math.imul(h, 0x846ca68b); h ^= h >>> 16;
        return ((h >>> 0) + 0.5) / 4294967296;
      };
      const gaussU = (a, b) => Math.sqrt(-2 * Math.log(hashU(a, 2 * b))) * Math.cos(2 * Math.PI * hashU(a, 2 * b + 1));
      let gas = null;
      function gasSetup() {
        const N = GAS_N, st = story, t0 = st.tShed;
        const g = { key: options.diskThick, t0, n: 0, next: 0, N,
          te: new Float64Array(N), origin: new Uint8Array(N), state: new Uint8Array(N),
          p: new Float64Array(3 * N), v: new Float64Array(3 * N), checkpoints: [], depositT: NaN };
        // Launch times: the inverse of the cumulative shedding rate.
        const M = 4000, cdf = new Float64Array(M + 1), span = st.tLaunchEnd - t0;
        for (let k = 0; k < M; k++) cdf[k + 1] = cdf[k] + st.rate(t0 + span * (k + 0.5) / M);
        let k = 0;
        for (let i = 0; i < N; i++) {
          const c = cdf[M] * (i + 0.5) / N;
          while (k < M - 1 && cdf[k + 1] < c) k++;
          g.te[i] = t0 + span * (k + (c - cdf[k]) / Math.max(cdf[k + 1] - cdf[k], 1e-300)) / M;
          g.origin[i] = hashU(i, 7777) < 0.6 ? 1 : 0;
        }
        g.checkpoints.push(gasSnapshot(g));
        return g;
      }
      function gasSnapshot(g) { return { n: g.n, next: g.next, p: Float32Array.from(g.p), v: Float32Array.from(g.v), state: g.state.slice() }; }
      function gasRestore(g, c) { g.n = c.n; g.next = c.next; g.p.set(c.p); g.v.set(c.v); g.state.set(c.state); }
      function gasLaunch(g, i, t) {
        const q = story.pairAt(t);
        if (!q) { g.state[i] = 2; return; }
        const k = g.origin[i], R = k ? bin.R2 : bin.R1;
        const xk = k ? q[2] : q[0], yk = k ? q[3] : q[1], xj = k ? q[0] : q[2], yj = k ? q[1] : q[3];
        const dx = xk - xj, dy = yk - yj, dl = Math.max(Math.hypot(dx, dy), 1e-6), ux = dx / dl, uy = dy / dl;
        const sp = 0.3 * R * gaussU(i, 1), sz = 0.25 * R * gaussU(i, 2);
        const x = xk + ux * 0.95 * R - uy * sp, y = yk + uy * 0.95 * R + ux * sp, z = sz;
        const w = story.rateAt(t), kick = 0.15 + 0.2 * hashU(i, 11);
        g.p[3 * i] = x; g.p[3 * i + 1] = y; g.p[3 * i + 2] = z;
        g.v[3 * i] = -w * y + ux * kick + 0.02 * gaussU(i, 3);
        g.v[3 * i + 1] = w * x + uy * kick + 0.02 * gaussU(i, 4);
        g.v[3 * i + 2] = 0.02 * gaussU(i, 5);
        g.state[i] = 1;
      }
      // One leapfrog step, from step n to n + 1.
      function gasStep(g) {
        const dt = GAS_DT, t = g.t0 + g.n * dt, t1 = t + dt, tm = bin.tMerge, Mt = bin.m1 + bin.m2;
        while (g.next < g.N && g.te[g.next] <= t) gasLaunch(g, g.next++, t);
        const qa = story.pairAt(t), qb = story.pairAt(t1);
        const R1 = 0.8 * bin.R1, R2 = 0.8 * bin.R2, rIn = 0.9 * remnantIsco;
        const post = t >= tm, Mc = qa ? Mt : bin.Mf, H = Math.max(options.diskThick, 0.06), n = g.n;
        const P = g.p, V = g.v;
        const acc = (q, x, y, z, out) => {
          if (q) {
            let dx = x - q[0], dy = y - q[1], r2 = dx * dx + dy * dy + z * z + 1.0, f1 = bin.m1 / (r2 * Math.sqrt(r2));
            let ex = x - q[2], ey = y - q[3], s2 = ex * ex + ey * ey + z * z + 1.0, f2 = bin.m2 / (s2 * Math.sqrt(s2));
            out[0] = -f1 * dx - f2 * ex; out[1] = -f1 * dy - f2 * ey; out[2] = -(f1 + f2) * z;
          } else {
            const r2 = x * x + y * y + z * z + 0.25, f = bin.Mf / (r2 * Math.sqrt(r2));
            out[0] = -f * x; out[1] = -f * y; out[2] = -f * z;
          }
        };
        const a = [0, 0, 0];
        for (let i = 0; i < g.next; i++) {
          if (g.state[i] !== 1) continue;
          const j = 3 * i;
          let x = P[j], y = P[j + 1], z = P[j + 2], vx = V[j], vy = V[j + 1], vz = V[j + 2];
          acc(qa, x, y, z, a); vx += 0.5 * dt * a[0]; vy += 0.5 * dt * a[1]; vz += 0.5 * dt * a[2];
          x += dt * vx; y += dt * vy; z += dt * vz;
          acc(qb, x, y, z, a); vx += 0.5 * dt * a[0]; vy += 0.5 * dt * a[1]; vz += 0.5 * dt * a[2];
          if (post) {
            const rc = Math.max(Math.hypot(x, y), 1e-3), om = Math.sqrt(Mc / (rc * rc * rc)), f = 1 - Math.exp(-dt * om / Math.PI);
            const ex = x / rc, ey = y / rc, vr = vx * ex + vy * ey, s = Math.sqrt(2 * f) * om * rc, seed = n * 3;
            vx += -vr * f * ex + 0.04 * s * gaussU(i + 100000, seed);
            vy += -vr * f * ey + 0.04 * s * gaussU(i + 100000, seed + 1);
            vz += -vz * f + H * s * gaussU(i + 100000, seed + 2);
          }
          let gone = x * x + y * y + z * z > 6400;
          if (qb) gone = gone || Math.hypot(x - qb[0], y - qb[1], z) < R1 || Math.hypot(x - qb[2], y - qb[3], z) < R2;
          else gone = gone || Math.hypot(x, y, z) < rIn;
          if (gone) { g.state[i] = 2; continue; }
          P[j] = x; P[j + 1] = y; P[j + 2] = z; V[j] = vx; V[j + 1] = vy; V[j + 2] = vz;
        }
        g.n++;
        if (g.n % GAS_EVERY === 0 && !g.checkpoints.some((c) => c.n === g.n)) g.checkpoints.push(gasSnapshot(g));
      }
      // The gas at time t, deposited into the texture.
      let gasData = null, gasTex = null, gasStatus = '';
      function gasAt(t) {
        if (!gas || gas.key !== options.diskThick) gas = gasSetup();
        const g = gas;
        if (g.depositT === t) return;
        const began = performance.now();
        const nT = Math.max(0, Math.floor((t - g.t0) / GAS_DT));
        if (nT < g.n) {
          let best = g.checkpoints[0];
          for (const c of g.checkpoints) if (c.n <= nT && c.n > best.n) best = c;
          gasRestore(g, best);
        }
        while (g.n < nT) gasStep(g);
        // Deposit, drifted from the step's time to t.
        const [nx, ny, nz] = GAS_GRID, out = options.diskOuter * unit;
        const B = 1.5 * Math.max(out, 12), Bz = 0.4 * B, vx = 2 * B / nx, vy = 2 * B / ny, vz = 2 * Bz / nz;
        if (!gasData) gasData = new Float32Array(4 * nx * ny * nz);
        gasData.fill(0);
        const dtD = t - (g.t0 + g.n * GAS_DT), fade = 0.8 * story.P;
        const sxy = 0.5 / (1.5 * 1.5), sz = 0.5 / (1.0 * 1.0);   // Gaussian widths in voxels: 1.5 across, 1 in z
        const wx = new Float64Array(5), wy = new Float64Array(5), wz = new Float64Array(3);
        for (let i = 0; i < g.next; i++) {
          if (g.state[i] !== 1) continue;
          const j = 3 * i;
          const fx = (g.p[j] + dtD * g.v[j] + B) / vx - 0.5, fy = (g.p[j + 1] + dtD * g.v[j + 1] + B) / vy - 0.5, fz = (g.p[j + 2] + dtD * g.v[j + 2] + Bz) / vz - 0.5;
          const ix = Math.round(fx), iy = Math.round(fy), iz = Math.round(fz);
          if (ix < 2 || iy < 2 || iz < 1 || ix > nx - 3 || iy > ny - 3 || iz > nz - 2) continue;
          let sx = 0, sy = 0, sw = 0;
          for (let k = 0; k < 5; k++) { const dx = ix + k - 2 - fx, dy = iy + k - 2 - fy; wx[k] = Math.exp(-dx * dx * sxy); wy[k] = Math.exp(-dy * dy * sxy); sx += wx[k]; sy += wy[k]; }
          for (let k = 0; k < 3; k++) { const dz = iz + k - 1 - fz; wz[k] = Math.exp(-dz * dz * sz); sw += wz[k]; }
          const norm = 1 / (sx * sy * sw);
          const young = 1 - smooth(t - g.te[i], 0, fade), y1 = g.origin[i] ? 0 : young, y2 = g.origin[i] ? young : 0;
          for (let c = 0; c < 3; c++) for (let b = 0; b < 5; b++) {
            const wbc = wz[c] * wy[b] * norm, row = 4 * (((iz + c - 1) * ny + (iy + b - 2)) * nx + (ix - 2));
            for (let a2 = 0; a2 < 5; a2++) {
              const wt = wbc * wx[a2], o = row + 4 * a2;
              gasData[o] += wt; gasData[o + 1] += wt * y1; gasData[o + 2] += wt * y2;
            }
          }
        }
        if (!gasTex) {
          gasTex = new THREE.Data3DTexture(gasData, nx, ny, nz);
          gasTex.format = THREE.RGBAFormat;
          gasTex.type = THREE.FloatType;
          gasTex.minFilter = gasTex.magFilter = renderer.extensions.has('OES_texture_float_linear') ? THREE.LinearFilter : THREE.NearestFilter;
          gasTex.unpackAlignment = 1;
        }
        gasTex.needsUpdate = true;
        uniforms.uGasTex.value = gasTex;
        uniforms.uGasMin.value.set(-B, -B, -Bz);
        uniforms.uGasSize.value.set(2 * B, 2 * B, 2 * Bz);
        // Optical depth per particle per unit volume.
        uniforms.uGasScale.value = GAS_KAPPA / (vx * vy * vz);
        g.depositT = t;
        let alive = 0, gone = 0;
        for (let i = 0; i < g.next; i++) { if (g.state[i] === 1) alive++; else if (g.state[i] === 2) gone++; }
        gasStatus = `gas: ${g.next} of ${g.N} particles shed, ${alive} in flight or orbit, ${gone} fallen back or accreted (${(performance.now() - began).toFixed(0)} ms)`;
      }
      // The disk's inner edge in units of the hole's mass (for its efficiency), and the last peak temperature.
      const rInOverM = bin ? remnantIsco / bin.Mf : M > 0 ? isco / M : 6;
      let diskPeak = 0, diskNormFor = -1;
      function diskReadout() {
        const nm = 2.8978e6 / diskPeak;   // Wien's displacement law
        const band = nm < 10 ? 'X-rays' : nm < 380 ? 'ultraviolet' : nm <= 750 ? 'visible' : 'infrared';
        const source = options.diskColour == 1 ? 'set directly'
          : `${ctx.fmt(Math.pow(10, options.diskMass))} M☉ at ${ctx.fmt(Math.pow(10, options.diskRate))} of the Eddington rate`;
        return `disk: peak ${ctx.fmt(diskPeak)} K (${source}); its spectrum peaks at ${nm >= 1000 ? ctx.fmt(nm / 1000) + ' µm' : ctx.fmt(nm) + ' nm'} (${band})`
          + (nm < 300 ? ', so visible light shows only the blue-white tail' : '');
      }

      const uniforms = {
        uCamPos: { value: new THREE.Vector3() }, uCamRight: { value: new THREE.Vector3() },
        uCamUp: { value: new THREE.Vector3() }, uCamFwd: { value: new THREE.Vector3() },
        uTanHalfFov: { value: 0.5 }, uAspect: { value: 1 },
        uRs: { value: params.rs }, uA: { value: params.a }, uQ2: { value: params.q * params.q },
        uHorizon: { value: horizon }, uEscape: { value: 100 }, uStepFrac: { value: options.accuracy },
        uMaxSteps: { value: options.steps }, uDisk: { value: 1 }, uShift: { value: 1 }, uSky: { value: 1 },
        uDiskIn: { value: isco || unit }, uDiskOut: { value: 8 * unit }, uJitter: { value: new THREE.Vector2() }, uFrameSeed: { value: 0 }, uPixelAngle: { value: 1e-3 },
        uDiskTPeak: { value: 5500 }, uDiskNorm: { value: 0 }, uDiskDetail: { value: 0.8 }, uDiskAge: { value: 0 },
        uDiskLayer: { value: [new THREE.Vector4(1, 0, 0, 0), new THREE.Vector4()] },
        uBB: { value: BB_TABLE.map((e) => new THREE.Vector4(e[0], e[1], e[2], e[3])) },
        uBodyCount: { value: 0 },
        uBody: { value: [new THREE.Vector4(), new THREE.Vector4(), new THREE.Vector4(), new THREE.Vector4()] },
        uBodyColor: { value: [new THREE.Color(), new THREE.Color(), new THREE.Color(), new THREE.Color()] },
        uLightMode: { value: 0 }, uSunDir: { value: new THREE.Vector3(0, 0, 1) }, uBodyGlow: { value: 0.06 }, uShadowR: { value: 1.5 * horizon },
        uMode: { value: 0 }, uTube: { value: null }, uTubeMin: { value: new THREE.Vector3() },
        uTubeSize: { value: new THREE.Vector3(1, 1, 1) }, uBodyReach: { value: 1e30 }, uMask: { value: null },
        uBinary: { value: 0 }, uX1: { value: new THREE.Vector3() }, uX2: { value: new THREE.Vector3() },
        uM1: { value: 0 }, uM2: { value: 0 }, uR1: { value: 0 }, uR2: { value: 0 }, uRrem: { value: 0 },
        uTides: { value: 0 }, uRemnant: { value: 0 }, uCM: { value: new THREE.Vector3() }, uL1: { value: new THREE.Vector3() },
        uOmega2: { value: 0 }, uLevel1: { value: -1e30 }, uLevel2: { value: -1e30 }, uB1: { value: 0 }, uB2: { value: 0 }, uGFloor: { value: 1e-3 },
        uBlend: { value: 0 }, uBall: { value: 0 }, uHeat: { value: 0 }, uBinW: { value: 1 }, uKerrW: { value: 0 },
        uDebris: { value: 0 }, uDebM: { value: 1 }, uDebA: { value: 0 }, uDebIn: { value: 1 }, uDebOut: { value: 10 }, uDebH: { value: 0.3 },
        uFill: { value: 0 }, uFront: { value: 0 },
        uGasTex: { value: null }, uGasMin: { value: new THREE.Vector3() }, uGasSize: { value: new THREE.Vector3(1, 1, 1) }, uGasScale: { value: 0 },
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
        uniforms.uFrameSeed.value = 97 * ((ctx.clock * 0.6180339887) % 1);
        // A binary: the first two tracks are the two masses (always traced; they are the spacetime), then the
        // remnant from the merger on. Its state is part of the camera key: everything is re-traced when it moves.
        let binKey = '';
        uniforms.uTides.value = 0; uniforms.uRemnant.value = 0; uniforms.uHeat.value = 0; uniforms.uBlend.value = 0; uniforms.uBall.value = 0;
        uniforms.uBinW.value = 1; uniforms.uKerrW.value = 0;
        uniforms.uDebris.value = 0; uniforms.uFill.value = 0;
        let gasMoves = false;
        if (bin) {
          const t = cLight * ctx.clock, merged = t >= bin.tMerge, tau = t - bin.tMerge, Mt = bin.m1 + bin.m2;
          if (story === undefined) story = buildStory();
          const stars = bin.R1 > 0 && bin.R2 > 0, tides = stars && !!options.tides;
          // The sketch of the merger itself (a remnant that collapses to a hole): see buildStory.
          const sketch = tides && bin.Rrem === 0 && !!(story && story.pairAt);
          let X1 = samplePhys(bodies[0].sample.x), X2 = samplePhys(bodies[1].sample.x);
          let binary = !merged, remnant = false, rs = merged ? 2 * bin.Mf : 0, spin = merged ? bin.aF : 0, horizonNow = merged ? remnantHorizon : 0;
          let shape = null, phase = '', kerrW = 0;
          if (tides && !merged) {
            shape = { omega2: story ? story.omega2At(t) : Mt / Math.pow(Math.max(X1.distanceTo(X2), 1e-6), 3), blend: 0, heat: -1 };
          } else if (sketch && tau < story.tEnd) {
            if (tau < story.tJoin) {
              // The pair closes up (story.pairAt).
              const u = tau / story.tJoin, q = story.pairAt(t);
              X1 = new THREE.Vector3(q[0], q[1], 0); X2 = new THREE.Vector3(q[2], q[3], 0);
              // The remnant's Kerr term fades in as the pair closes up (see binaryMetric).
              binary = true; rs = 2 * bin.Mf; spin = bin.aF; horizonNow = 0;
              kerrW = smooth(u, 0.05, 1);
              shape = { omega2: story.omega * story.omega, blend: smooth(u, 0.5, 0.95), heat: story.heatM + (0.8 - story.heatM) * smooth(u, 0, 1), c: 1 };
              phase = `merged: the remnant spins and closes into one body (${Math.round(100 * u)}%)`;
            } else {
              // The remnant's Kerr metric; the ball shrinks onto the horizon.
              const v = (tau - story.tJoin) / (story.tEnd - story.tJoin);
              binary = false; remnant = true;
              X1 = new THREE.Vector3(); X2 = new THREE.Vector3();
              shape = { blend: 1, heat: 0.8, ball: Math.max(story.rJoin + (1.02 * horizonNow - story.rJoin) * smooth(v, 0.05, 0.95), 1.02 * horizonNow) };
              phase = `merged: the remnant collapses into the hole (${Math.round(100 * v)}%)`;
            }
          }
          uniforms.uBinary.value = binary ? 1 : 0;
          uniforms.uRemnant.value = remnant ? 1 : 0;
          uniforms.uKerrW.value = kerrW; uniforms.uBinW.value = 1 - kerrW;
          uniforms.uM1.value = bin.m1; uniforms.uM2.value = bin.m2;
          uniforms.uR1.value = bin.R1; uniforms.uR2.value = bin.R2;
          uniforms.uX1.value.copy(X1); uniforms.uX2.value.copy(X2);
          const star = new THREE.Color(0.85, 0.92, 1.0);
          uniforms.uStarColour1.value.copy(star).lerp(bodies[0].color, 0.3).multiplyScalar(1.5);
          uniforms.uStarColour2.value.copy(star).lerp(bodies[1].color, 0.3).multiplyScalar(1.5);
          uniforms.uRs.value = rs;
          uniforms.uA.value = spin;
          uniforms.uHorizon.value = binary ? 0 : horizonNow;
          uniforms.uRrem.value = merged && !shape ? bin.Rrem : 0;
          uniforms.uDiskIn.value = merged ? remnantIsco : unit;
          // No disk about the pair's centre of mass: a disk only around the remnant (with the sketch, its torus).
          if (!merged || sketch) uniforms.uDisk.value = 0;
          const sep = X1.distanceTo(X2);
          binaryState = phase || (merged ? `merged: ${bin.Rrem > 0 ? 'star' : 'Kerr hole'} M = ${ctx.fmt(bin.Mf)}, a/M = ${ctx.fmt(bin.aF / bin.Mf)}` : `inspiral: separation ${ctx.fmt(sep)}`);
          if (shape) {
            uniforms.uTides.value = 1;
            uniforms.uBlend.value = shape.blend;
            if (shape.blend < 0.999) {
              const rl = roche(X1, X2, bin.m1, bin.m2, bin.R1, bin.R2, shape.omega2, shape.c !== undefined ? { c: shape.c } : {});
              uniforms.uCM.value.copy(rl.cm); uniforms.uL1.value.copy(rl.L1); uniforms.uOmega2.value = shape.omega2;
              uniforms.uLevel1.value = rl.level1; uniforms.uLevel2.value = rl.level2; uniforms.uB1.value = rl.B1; uniforms.uB2.value = rl.B2;
              uniforms.uGFloor.value = 0.3 * Math.min(bin.m1 / (bin.R1 * bin.R1), bin.m2 / (bin.R2 * bin.R2));
              uniforms.uBall.value = rl.rEq;
              if (shape.heat < 0) {
                // Before the merger: matter at the neck heats as the surfaces close on the inner Lagrange point.
                shape.heat = 0.35 * smooth(rl.q, -0.05, 0.15);
                if (rl.connected) binaryState += ' (in contact)';
                else if (story && story.tShed !== undefined && t >= story.tShed) binaryState += ' (shedding matter from the outer tips)';
              }
            } else uniforms.uBall.value = shape.ball ?? story.rJoin;
            uniforms.uHeat.value = shape.heat;
          }
          if (sketch && options.disk && t >= story.tShed) {
            // The gas: the particles shed by the stars, at this time (gasAt), as a density texture.
            uniforms.uDebris.value = 1;
            uniforms.uDebM.value = bin.Mf; uniforms.uDebA.value = bin.aF; uniforms.uDebIn.value = remnantIsco;
            uniforms.uDebH.value = Math.max(options.diskThick, 0.06); uniforms.uDebOut.value = options.diskOuter * unit;
            gasAt(t);
            gasMoves = merged;
            if (merged && !phase && tau < 5 * story.P) binaryState += ' (the debris settles into a torus)';
          }
          const timed = !!shape || uniforms.uDebris.value > 0;
          binKey = [merged, X1.toArray().map(v => v.toPrecision(6)).join(','), X2.toArray().map(v => v.toPrecision(6)).join(','), timed ? t.toPrecision(7) : ''].join('/');
        }
        // A thick disk otherwise (option): the same gas volume, filled, about the hole.
        if (uniforms.uDebris.value === 0 && uniforms.uDisk.value === 1 && options.diskThick > 0) {
          uniforms.uDebris.value = 2;
          uniforms.uDisk.value = 0;
          uniforms.uDebM.value = 0.5 * uniforms.uRs.value; uniforms.uDebA.value = uniforms.uA.value;
          uniforms.uDebIn.value = uniforms.uDiskIn.value; uniforms.uDebOut.value = uniforms.uDiskOut.value; uniforms.uDebH.value = options.diskThick;
          uniforms.uFill.value = 1; uniforms.uFront.value = 2 * uniforms.uDiskOut.value;
        }
        // The disk: peak temperature and its brightness unit; the pattern's shear age (4 radians of turn at the
        // inner edge) and, when it rotates, two layers that each live one inner orbit and cross-fade, so the
        // winding stays bounded while every ring turns at its own speed. A merger's debris always turns.
        diskPeak = options.diskColour == 1 ? options.diskTemp : diskPeakTemperature(Math.pow(10, options.diskMass), Math.pow(10, options.diskRate), rInOverM);
        if (diskPeak !== diskNormFor) { uniforms.uDiskNorm.value = Math.log10(diskMeanLuminance(diskPeak)); diskNormFor = diskPeak; }
        uniforms.uDiskTPeak.value = diskPeak;
        uniforms.uDiskDetail.value = options.diskDetail;
        const gasOn = uniforms.uDebris.value > 0;
        const Md = gasOn ? uniforms.uDebM.value : 0.5 * uniforms.uRs.value, rIn = gasOn ? uniforms.uDebIn.value : uniforms.uDiskIn.value;
        const spinIn = gasOn ? uniforms.uDebA.value : uniforms.uA.value, Mr = Math.max(Md * rIn - (gasOn ? 0 : uniforms.uQ2.value), 0);
        const omegaIn = Mr > 0 ? Math.sqrt(Mr) / (rIn * rIn + spinIn * Math.sqrt(Mr)) : 0;
        const ageFull = omegaIn > 0 ? 4 / omegaIn : 0;
        uniforms.uDiskAge.value = ageFull;
        const layers = uniforms.uDiskLayer.value;
        let diskKey = '';
        if (((options.disk && options.diskMotion) || gasMoves) && omegaIn > 0) {
          // (A merger's torus: two inner orbits per layer, so it cross-fades less between the frames of a video.)
          const period = (gasMoves ? 2 : 1) * 2 * Math.PI / omegaIn, t = cLight * ctx.clock;
          for (let k = 0; k < 2; k++) {
            const u = t / period + k / 2, s = u - Math.floor(u);
            layers[k].set(Math.sin(Math.PI * s) ** 2, (s - 0.5) * period, Math.floor(u) + 13 * k, 0);
          }
          diskKey = t.toPrecision(7);
        } else {
          layers[0].set(1, 0, 0, 0);
          layers[1].set(0, 0, 0, 0);
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
        {
          uniforms.uLightMode.value = Math.round(options.bodyLight);
          const az = THREE.MathUtils.degToRad(options.sunAz), el = THREE.MathUtils.degToRad(options.sunEl);
          uniforms.uSunDir.value.set(Math.cos(el) * Math.cos(az), Math.cos(el) * Math.sin(az), Math.sin(el));
          uniforms.uBodyGlow.value = options.bodyGlow;
          uniforms.uShadowR.value = 1.5 * uniforms.uHorizon.value;
        }
        // What changes the image: the camera key (camera, options, window, which bodies are shown) restarts
        // everything; the body key (where they are) only the body layer. The camera is rounded so that motion
        // far below a pixel (the tail of an easing camera) does not restart refinement.
        const { exposure, paths, ...traced } = options;
        return {
          cam: [e.map(v => v.toPrecision(5)).join(','), shownIndex.join(','), JSON.stringify(traced), renderer.domElement.width, renderer.domElement.height, binKey, diskKey].join('|'),
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
        uniforms.uPixelAngle.value = 2 * uniforms.uTanHalfFov.value / targetRT.height;
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

      // A merger scene starts at the chosen separation (option startSep).
      function applyStart() {
        if (!bin) return;
        if (story === undefined) story = buildStory();
        if (!story) return;
        ctx.setTimeRange(story.startAt(options.startSep) / cLight, ctx.tEnd);
      }
      applyStart();

      return {
        group: new THREE.Group(),
        hidesBodies: true,
        // For the Render dialog: the span around a merger that it can play in slow motion (clock units).
        slowWindow() { return story && story.slow ? { t0: story.slow[0] / cLight, t1: story.slow[1] / cLight, label: 'the merger', factor: story.slowFactor } : null; },
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
          if (id === 'startSep') { applyStart(); return true; }
          return true;   // every other option is a uniform, picked up on the next frame
        },
        // For offline rendering: the frame is final (all samples, including the body layer's).
        converged() {
          if (level !== 3) return false;
          if (floatOk && bodyCount > 0 && maskReady) return dynFull && dynSample >= Math.max(1, Math.round(options.samples));
          return true;
        },
        readout() {
          return [bin ? binaryState : `horizon r+ ${ctx.fmt(horizon)}, ISCO ${ctx.fmt(isco)}`, status]
            .concat(options.disk ? [diskReadout()] : [], uniforms.uDebris.value === 1 && gasStatus ? [gasStatus] : [], dynStatus ? [dynStatus] : []);
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
          if (gasTex) gasTex.dispose();
          material.dispose();
          addMaterial.dispose();
          displayMaterial.dispose();
        },
      };
    },
  });
})();
