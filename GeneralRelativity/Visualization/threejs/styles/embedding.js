// Embedding style: the equatorial plane as a surface of revolution (Flamm's paraboloid for
// Schwarzschild) coloured by the Kretschmann scalar, reference spheres, scalar field slices, and the
// bodies in 3D. Works for every scene; without a surface it is a plain 3D view.
(window.GRStyles = window.GRStyles || []).push({
  id: 'embedding',
  label: 'Embedding / 3D',
  order: 1,
  description: 'The equatorial plane embedded in Euclidean space (depth shows the stretching of radial distance), coloured by curvature; spheres and field slices in 3D.',
  options: [
    { id: 'surface', label: 'surface', type: 'checkbox', value: 1 },
    { id: 'wire', label: 'wireframe', type: 'checkbox', value: 1 },
    { id: 'fields', label: 'field slices', type: 'checkbox', value: 1 },
  ],
  create(ctx) {
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
