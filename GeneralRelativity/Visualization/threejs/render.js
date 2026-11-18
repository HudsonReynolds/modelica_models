// Offline rendering: steps the simulation clock through a time range with the camera fixed, waits until the
// active style reports each frame as final (the ray tracer: every refinement sample traced), and writes the
// frames as a video or an image sequence. A style with a key event (style.slowWindow(), a merger) can be
// rendered in slow motion around it: frames denser in simulation time by a factor there.
//
// Outputs:
//   - in the browser (any page, also opened from disk): MP4 (H.264) or WebM (VP9) encoded with WebCodecs and
//     muxed by mp4-muxer / webm-muxer, or a .zip of PNG frames; the file is saved through the browser's
//     download (a save dialog where the browser asks where to save, or showSaveFilePicker where supported);
//   - through the local helper render_server.py (when the viewer is opened from it, http://127.0.0.1:8765):
//     PNG frames written to a folder typed in the dialog, then MPEG-2 (.mpeg) or H.264 (.mp4) made by ffmpeg.
//
// Classic script (not a module) so the viewer still opens from disk; the core calls GRRender.install(api).
(function () {
  const MP4_MUXER = 'https://cdn.jsdelivr.net/npm/mp4-muxer@5.2.2/+esm';
  const WEBM_MUXER = 'https://cdn.jsdelivr.net/npm/webm-muxer@5.1.4/+esm';
  const RESOLUTIONS = [[1280, 720], [1920, 1080], [2560, 1440], [3840, 2160], [1080, 1080], [1080, 1920]];
  // Frame rates MPEG-2 allows.
  const MPEG_RATES = [23.976, 24, 25, 29.97, 30, 50, 59.94, 60];

  const css = `
    #renderDialog { background: #151a24; color: #e6e9ef; border: 1px solid #2f3a4e; border-radius: 8px; padding: 14px 16px; width: min(560px, calc(100vw - 32px)); }
    #renderDialog::backdrop { background: rgba(0, 0, 0, 0.45); }
    #renderDialog h3 { margin: 0 0 8px; font-size: 15px; }
    #renderDialog .grid { display: grid; grid-template-columns: 150px minmax(0, 1fr); gap: 6px 10px; align-items: center; font-size: 13px; }
    #renderDialog .grid > label { color: #8a93a6; }
    #renderDialog input[type=number], #renderDialog input[type=text] { background: #0f131b; color: #e6e9ef; border: 1px solid #2f3a4e; border-radius: 4px; padding: 3px 6px; font: inherit; width: 100%; box-sizing: border-box; }
    #renderDialog .row { display: flex; gap: 6px; align-items: center; }
    #renderDialog .row input[type=number] { width: 90px; }
    #renderDialog .note { color: #8a93a6; font-size: 12px; margin: 8px 0 0; white-space: pre-line; }
    #renderDialog .warn { color: #f0c674; }
    #renderDialog .buttons { display: flex; gap: 8px; justify-content: flex-end; margin-top: 12px; }
    #renderDialog progress { width: 100%; height: 10px; margin-top: 10px; }
    #renderDialog .status { font-family: ui-monospace, monospace; font-size: 12px; white-space: pre-line; margin-top: 6px; }
    #renderDialog .hidden { display: none; }
    #renderDialog input.command { font-family: ui-monospace, monospace; font-size: 11px; }
  `;

  // --- small helpers ---------------------------------------------------------------------------------
  const el = (tag, props = {}, children = []) => {
    const e = document.createElement(tag);
    Object.assign(e, props);
    for (const c of [].concat(children)) e.append(c);
    return e;
  };
  const even = (n) => Math.max(2, 2 * Math.round(n / 2));
  const download = (blob, name) => {
    const a = el('a', { href: URL.createObjectURL(blob), download: name });
    document.body.append(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(a.href), 120000);
  };
  const clockText = (s) => {
    if (!Number.isFinite(s)) return '?';
    const m = Math.floor(s / 60), r = Math.round(s - 60 * m);
    return m > 0 ? `${m} min ${r} s` : `${r} s`;
  };
  const safeName = (s) => String(s).replace(/[^A-Za-z0-9._-]+/g, '_').replace(/_+/g, '_').replace(/^_+|_+$/g, '') || 'render';

  // CRC-32 (for the stored, uncompressed ZIP of PNG frames).
  const CRC_TABLE = (() => {
    const t = new Uint32Array(256);
    for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; }
    return t;
  })();
  const crc32 = (bytes) => {
    let c = 0xffffffff;
    for (let i = 0; i < bytes.length; i++) c = CRC_TABLE[(c ^ bytes[i]) & 0xff] ^ (c >>> 8);
    return (c ^ 0xffffffff) >>> 0;
  };

  // Detect the helper server (only when the page itself comes from it).
  async function helperInfo() {
    if (!/^https?:$/.test(location.protocol)) return null;
    try {
      const r = await fetch('/api/info', { cache: 'no-store' });
      if (!r.ok) return null;
      const info = await r.json();
      return info && info.token ? info : null;
    } catch (e) { return null; }
  }

  async function h264Codec(width, height, fps, bitrate) {
    const px = width * height;
    // Profile High; level by frame size, with fallbacks.
    const levels = px <= 921600 ? ['64001f', '640028', '640033'] : px <= 2088960 ? ['640028', '64002a', '640033'] : px <= 3686400 ? ['640032', '640033'] : ['640033', '640034'];
    for (const l of levels.concat(['4d0033', '42e033'])) {
      const config = { codec: 'avc1.' + l, width, height, bitrate, framerate: fps, avc: { format: 'avc' } };
      try { if ((await VideoEncoder.isConfigSupported(config)).supported) return config; } catch (e) { /* try the next */ }
    }
    return null;
  }
  async function vp9Codec(width, height, fps, bitrate) {
    for (const codec of ['vp09.00.41.08', 'vp09.00.10.08', 'vp8']) {
      const config = { codec, width, height, bitrate, framerate: fps };
      try { if ((await VideoEncoder.isConfigSupported(config)).supported) return config; } catch (e) { /* next */ }
    }
    return null;
  }

  // --- output sinks: capture(canvas, k) runs synchronously right after the final draw of frame k; flush()
  // finishes that frame (encoding back-pressure, PNG upload); finish() writes the result. -----------------
  async function encoderSink(cfg, kind) {
    const bitrate = Math.round(cfg.mbps * 1e6);
    const config = kind === 'mp4' ? await h264Codec(cfg.width, cfg.height, cfg.fps, bitrate) : await vp9Codec(cfg.width, cfg.height, cfg.fps, bitrate);
    if (!config) throw new Error(`This browser cannot encode ${kind === 'mp4' ? 'H.264' : 'VP9'} at ${cfg.width} x ${cfg.height}.`);
    const lib = await import(kind === 'mp4' ? MP4_MUXER : WEBM_MUXER);
    let fileStream = null;
    if (cfg.fileHandle) fileStream = await cfg.fileHandle.createWritable();
    const target = fileStream ? new lib.FileSystemWritableFileStreamTarget(fileStream) : new lib.ArrayBufferTarget();
    const video = kind === 'mp4'
      ? { codec: 'avc', width: cfg.width, height: cfg.height, frameRate: cfg.fps }
      : { codec: config.codec === 'vp8' ? 'V_VP8' : 'V_VP9', width: cfg.width, height: cfg.height, frameRate: cfg.fps };
    const muxer = new lib.Muxer(kind === 'mp4' ? { target, video, fastStart: fileStream ? false : 'in-memory' } : { target, video });
    let failure = null;
    const encoder = new VideoEncoder({ output: (chunk, meta) => muxer.addVideoChunk(chunk, meta), error: (e) => { failure = e; } });
    encoder.configure(config);
    const us = 1e6 / cfg.fps;
    return {
      describe: `${kind === 'mp4' ? 'MP4' : 'WebM'} (${config.codec}, ${cfg.mbps} Mbit/s)`,
      capture(canvas, k) {
        const frame = new VideoFrame(canvas, { timestamp: Math.round(k * us), duration: Math.round(us) });
        encoder.encode(frame, { keyFrame: k % Math.max(1, Math.round(2 * cfg.fps)) === 0 });
        frame.close();
      },
      async flush() {
        while (encoder.encodeQueueSize > 2) await new Promise((r) => setTimeout(r, 4));
        if (failure) throw failure;
      },
      async finish() {
        await encoder.flush();
        if (failure) throw failure;
        muxer.finalize();
        encoder.close();
        if (fileStream) { await fileStream.close(); return `saved to ${cfg.fileHandle.name}`; }
        const blob = new Blob([target.buffer], { type: kind === 'mp4' ? 'video/mp4' : 'video/webm' });
        download(blob, cfg.fileName);
        return `${(blob.size / 1e6).toFixed(1)} MB, saved through the browser as ${cfg.fileName}`;
      },
      abort() { try { encoder.close(); } catch (e) { /* closed */ } if (fileStream) fileStream.abort().catch(() => {}); },
    };
  }

  function zipSink(cfg) {
    // A stored (uncompressed) ZIP: PNG data is already compressed. Parts are kept as Blobs, which browsers may
    // keep on disk, so long sequences do not need one contiguous buffer.
    const parts = [], central = [];
    let offset = 0, pending = null;
    const u16 = (v) => [v & 0xff, (v >>> 8) & 0xff];
    const u32 = (v) => [v & 0xff, (v >>> 8) & 0xff, (v >>> 16) & 0xff, (v >>> 24) & 0xff];
    const base = cfg.fileName.replace(/\.zip$/i, '');
    return {
      describe: 'PNG frames in a .zip',
      capture(canvas, k) { pending = new Promise((res) => canvas.toBlob(res, 'image/png')).then((b) => ({ blob: b, k })); },
      async flush() {
        const { blob, k } = await pending;
        const bytes = new Uint8Array(await blob.arrayBuffer());
        const name = new TextEncoder().encode(`${base}/frame_${String(k).padStart(5, '0')}.png`);
        const crc = crc32(bytes);
        const header = new Uint8Array([0x50, 0x4b, 0x03, 0x04, ...u16(20), ...u16(0), ...u16(0), ...u16(0), ...u16(0x21),
          ...u32(crc), ...u32(bytes.length), ...u32(bytes.length), ...u16(name.length), ...u16(0)]);
        parts.push(new Blob([header, name, bytes]));
        central.push({ name, crc, size: bytes.length, offset });
        offset += header.length + name.length + bytes.length;
      },
      async finish() {
        const dir = [];
        let size = 0;
        for (const c of central) {
          const h = new Uint8Array([0x50, 0x4b, 0x01, 0x02, ...u16(20), ...u16(20), ...u16(0), ...u16(0), ...u16(0), ...u16(0x21),
            ...u32(c.crc), ...u32(c.size), ...u32(c.size), ...u16(c.name.length), ...u16(0), ...u16(0), ...u16(0), ...u16(0), ...u32(0), ...u32(c.offset)]);
          dir.push(h, c.name);
          size += h.length + c.name.length;
        }
        const end = new Uint8Array([0x50, 0x4b, 0x05, 0x06, ...u16(0), ...u16(0), ...u16(central.length), ...u16(central.length), ...u32(size), ...u32(offset), ...u16(0)]);
        const blob = new Blob([...parts, ...dir, end], { type: 'application/zip' });
        download(blob, cfg.fileName);
        return `${central.length} frames, ${(blob.size / 1e6).toFixed(1)} MB, saved through the browser as ${cfg.fileName}`;
      },
      abort() {},
    };
  }

  async function helperSink(cfg, helper) {
    const post = async (path, body, type = 'application/json') => {
      const r = await fetch(path, { method: 'POST', headers: { 'Content-Type': type, 'X-Render-Token': helper.token }, body });
      const reply = await r.json().catch(() => ({ ok: false, error: `HTTP ${r.status}` }));
      if (!r.ok || !reply.ok) throw new Error(reply.error || `HTTP ${r.status}`);
      return reply;
    };
    const start = await post('/api/render/start', JSON.stringify({ folder: cfg.folder, name: cfg.baseName, format: cfg.format, fps: cfg.fps, width: cfg.width, height: cfg.height, keepFrames: cfg.keepFrames }));
    let pending = null;
    return {
      describe: { 'h-mpeg': 'MPEG-2 .mpeg (ffmpeg)', 'h-mp4': 'H.264 .mp4 (ffmpeg)', 'h-png': 'PNG frames' }[cfg.format] + ` in ${start.folder}`,
      capture(canvas, k) { pending = new Promise((res) => canvas.toBlob(res, 'image/png')).then((b) => ({ blob: b, k })); },
      async flush() { const { blob, k } = await pending; await post(`/api/render/frame?i=${k}`, blob, 'image/png'); },
      async finish() { const done = await post('/api/render/finish', '{}'); return done.message; },
      abort() { post('/api/render/cancel', '{}').catch(() => {}); },
    };
  }

  // --- the dialog and the render loop ---------------------------------------------------------------
  window.GRRender = {
    install(api) {
      document.head.append(el('style', { textContent: css }));
      const button = el('button', { textContent: 'Render…', title: 'Render the animation to a video or image sequence, frame by frame with the camera fixed' });
      api.header.append(button);
      const dialog = el('dialog', { id: 'renderDialog' });
      document.body.append(dialog);
      let helper = null, job = null;

      const field = (label, input, help) => [el('label', { textContent: label, title: help || '' }), input];
      const number = (value, step = 'any', min = '') => el('input', { type: 'number', value, step, min });
      const resolution = el('select');
      const width = number(1920, 2, 16), height = number(1080, 2, 16);
      const fps = number(24, 'any', 1);
      const t0 = number(0), t1 = number(1);
      const seconds = number(10, 'any', 0.05);
      const frameInfo = el('span', { className: 'note' });
      const samples = number(8, 1, 1);
      const mbps = number(12, 'any', 0.5);
      const format = el('select');
      const fileName = el('input', { type: 'text' });
      const chooseFile = el('button', { textContent: 'Choose…', type: 'button', title: 'Pick where to save (browsers with showSaveFilePicker)' });
      const folder = el('input', { type: 'text' });
      const keepFrames = el('input', { type: 'checkbox' });
      const diskMotion = el('input', { type: 'checkbox' });
      const slow = number(4, 'any', 1);
      const destNote = el('div', { className: 'note' });
      const command = el('input', { type: 'text', readOnly: true, className: 'command' });
      const copyCommand = el('button', { textContent: 'Copy', type: 'button', title: 'Copy the command that renders this from a terminal, without the browser (render_cli.py)' });
      const progress = el('progress', { max: 1, value: 0, className: 'hidden' });
      const status = el('div', { className: 'status' });
      const start = el('button', { textContent: 'Render', type: 'button' });
      const close = el('button', { textContent: 'Close', type: 'button' });
      let fileHandle = null;

      const rows = {
        size: el('div', { className: 'row' }, [width, '×', height]),
        time: el('div', { className: 'row' }, [t0, 'to', t1]),
        file: el('div', { className: 'row' }, [fileName, chooseFile]),
        keep: el('div', { className: 'row' }, [keepFrames, el('span', { textContent: 'also keep the PNG frames' })]),
        disk: el('div', { className: 'row' }, [diskMotion, el('span', { textContent: 'the accretion disk rotates' })]),
        slow: el('div', { className: 'row' }, [slow, el('span', { className: 'slowLabel', textContent: '× slower around the merger' })]),
        command: el('div', { className: 'row' }, [command, copyCommand]),
      };
      const grid = el('div', { className: 'grid' }, [
        ...field('resolution', resolution), ...field('', rows.size),
        ...field('frame rate (fps)', fps),
        ...field('simulation time', rows.time, 'Clock range to render (the scene clock, as in the side panel)'),
        ...field('video length (s)', seconds),
        ...field('slow motion', rows.slow, 'Frames are this many times denser in simulation time around the key event (a neutron-star merger: from the first shedding through the collapse and the torus forming; black holes: the last orbits and the ringdown), easing in and out; 1 keeps the clock steady. Renders only, not live playback'),
        ...field('', frameInfo),
        ...field('samples per frame', samples, 'Ray tracer: jittered samples averaged per frame (anti-aliasing). 8 is indistinguishable from 48 in a video (mean difference 0.2 of 255 levels, far below video compression); more for stills'),
        ...field('', rows.disk, 'Ray tracer: the disk pattern turns at the orbital speed of each radius (the style option "disk rotates", for this render)'),
        ...field('format', format),
        ...field('bit rate (Mbit/s)', mbps, 'Browser encoders (MP4, WebM)'),
        ...field('file name', rows.file),
        ...field('output folder', folder, 'Folder on this computer (helper server); ~ is your home folder'),
        ...field('', rows.keep),
        ...field('command line', rows.command, 'The same render from a terminal, in a headless browser (python3 render_cli.py, in the viewer folder): it keeps going with the browser closed, and can be split across machines'),
      ]);
      dialog.append(el('h3', { textContent: 'Render animation' }), grid, destNote, progress, status, el('div', { className: 'buttons' }, [start, close]));
      const rowOf = (input) => input.closest('.row') || input;
      const hasDisk = () => { const o = api.options(); return !!(o && 'diskMotion' in o && o.disk); };
      const show = (input, on) => { const node = rowOf(input); node.classList.toggle('hidden', !on); node.previousElementSibling.classList.toggle('hidden', !on); };

      // The command that renders the same thing from a terminal (render_cli.py): scene, style, the options that
      // differ from their defaults, bodies hidden, legend, camera, time range, frames and output.
      const shellQuote = (v) => (/^[A-Za-z0-9._\/:=,+-]+$/.test(String(v)) ? String(v) : `'${String(v).replace(/'/g, `'\\''`)}'`);
      const num = (v) => String(+Number(v).toPrecision(7));
      function commandLine() {
        // Run in the viewer folder (where render_cli.py and, normally, the scenes are), when the page knows it.
        const dir0 = helper && helper.root ? helper.root : location.protocol === 'file:' ? decodeURIComponent(location.pathname.replace(/\/[^/]*$/, '')) : '';
        const parts = (dir0 ? ['cd', shellQuote(dir0), '&&'] : []).concat(['python3', 'render_cli.py', shellQuote(api.sceneFile() || 'SCENE.ndjson'), '--style', api.styleId()]);
        const o = api.options() || {}, d = api.optionDefaults();
        for (const [k, v] of Object.entries(o)) if (k !== 'samples' && k !== 'diskMotion' && d[k] !== undefined && v !== d[k]) parts.push('--option', shellQuote(`${k}=${v}`));
        for (const b of api.hiddenBodies()) parts.push('--hide', shellQuote(b));
        if (!api.legendShown()) parts.push('--no-legend');
        const c = api.camera.position, tg = api.controls.target;
        parts.push('--camera', [c.x, c.y, c.z].map(num).join(','), '--target', [tg.x, tg.y, tg.z].map(num).join(','));
        parts.push('--from', num(t0.value), '--to', num(t1.value), '--fps', num(fps.value), '--seconds', num(seconds.value));
        parts.push('--size', `${even(Number(width.value))}x${even(Number(height.value))}`);
        if (api.options() && 'samples' in api.options()) parts.push('--samples', String(Math.max(1, Math.round(Number(samples.value)))));
        if (slowWindow()) parts.push('--slow', num(slow.value));
        if (hasDisk() && diskMotion.checked) parts.push('--disk-rotates');
        const f = format.value, base = safeName(fileName.value.replace(/\.(mp4|webm|zip|mpeg)$/i, ''));
        const ext = f === 'h-mpeg' ? '.mpeg' : f === 'h-png' || f === 'zip' ? '/' : '.mp4';
        const dir = (folder.value || '~/Documents/Visualization/renders').replace(/\/+$/, '');
        parts.push('-o', (dir.startsWith('~/') ? '~/' + shellQuote(dir.slice(2)) : shellQuote(dir)) + '/' + base + ext);
        if (keepFrames.checked && f !== 'h-png' && f !== 'zip') parts.push('--keep-frames');
        return parts.join(' ');
      }
      copyCommand.addEventListener('click', async () => {
        command.value = commandLine();
        command.select();
        try { await navigator.clipboard.writeText(command.value); copyCommand.textContent = 'Copied'; }
        catch (e) { document.execCommand && document.execCommand('copy'); copyCommand.textContent = 'Copied'; }
        setTimeout(() => { copyCommand.textContent = 'Copy'; }, 1500);
      });

      const formats = () => {
        const list = [];
        if (helper) {
          list.push({ value: 'h-mpeg', label: 'MPEG-2 .mpeg (ffmpeg, to folder)' }, { value: 'h-mp4', label: 'H.264 .mp4 (ffmpeg, to folder)' }, { value: 'h-png', label: 'PNG frames (to folder)' });
        }
        if (typeof VideoEncoder === 'function') list.push({ value: 'mp4', label: 'MP4, H.264 (browser)' }, { value: 'webm', label: 'WebM, VP9 (browser)' });
        list.push({ value: 'zip', label: 'PNG frames, .zip (browser)' });
        return list;
      };
      const extension = (f) => ({ mp4: '.mp4', webm: '.webm', zip: '.zip', 'h-mpeg': '.mpeg', 'h-mp4': '.mp4', 'h-png': '' }[f]);
      const frames = () => Math.max(1, Math.round(Number(fps.value) * Number(seconds.value)));
      // Slow motion: the style's window [t0, t1] around its key event (style.slowWindow()), where frames are
      // S times denser in simulation time, with smooth ramps of half the window's length on either side. The
      // frame times invert the cumulative frame density; S = 1 (or no window) gives evenly spaced frames.
      const slowWindow = () => { const st = api.style(); return st && st.slowWindow ? st.slowWindow() : null; };
      function frameTimes(a, b, n, S) {
        const lin = (k) => (n > 1 ? a + (b - a) * k / (n - 1) : a);
        const win = slowWindow();
        if (!(S > 1) || !win || n < 2 || !(b > a)) return Array.from({ length: n }, (_, k) => lin(k));
        const ramp = 0.5 * (win.t1 - win.t0);
        const ss = (x, e0, e1) => { const t = Math.min(Math.max((x - e0) / (e1 - e0), 0), 1); return t * t * (3 - 2 * t); };
        const density = (t) => 1 + (S - 1) * ss(t, win.t0 - ramp, win.t0) * (1 - ss(t, win.t1, win.t1 + ramp));
        const M = 8192, h = (b - a) / M, C = new Float64Array(M + 1);
        for (let i = 0; i < M; i++) C[i + 1] = C[i] + 0.5 * h * (density(a + i * h) + density(a + (i + 1) * h));
        const out = [];
        let i = 0;
        for (let k = 0; k < n; k++) {
          const c = C[M] * k / (n - 1);
          while (i < M - 1 && C[i + 1] < c) i++;
          out.push(k === n - 1 ? b : a + h * (i + (c - C[i]) / Math.max(C[i + 1] - C[i], 1e-300)));
        }
        return out;
      }
      const refresh = () => {
        const custom = resolution.value === 'custom';
        if (!custom) { const [w, h] = resolution.value.split('x').map(Number); width.value = w; height.value = h; }
        width.disabled = height.disabled = !custom;
        const n = frames(), span = Number(t1.value) - Number(t0.value);
        const win = slowWindow(), S = Number(slow.value);
        show(slow, !!win);
        if (win) rows.slow.querySelector('.slowLabel').textContent = `× slower around ${win.label || 'the key event'}`;
        if (win && S > 1 && n > 2) {
          const ts = frameTimes(Number(t0.value), Number(t1.value), n, S);
          let lo = Infinity, hi = 0;
          for (let k = 1; k < n; k++) { const d = ts[k] - ts[k - 1]; lo = Math.min(lo, d); hi = Math.max(hi, d); }
          frameInfo.textContent = `${n} frames; ${api.fmt(hi)} of simulation time per frame, down to ${api.fmt(lo)} around ${win.label || 'the key event'}`;
        } else frameInfo.textContent = `${n} frames; ${api.fmt(n > 1 ? span / (n - 1) : 0)} of simulation time per frame, ${api.fmt(span / Math.max(Number(seconds.value), 1e-9))} per second of video`;
        const f = format.value, viaHelper = f.startsWith('h-');
        show(folder, viaHelper);
        show(keepFrames, viaHelper && f !== 'h-png');
        show(mbps, f === 'mp4' || f === 'webm');
        chooseFile.classList.toggle('hidden', !(window.showSaveFilePicker && (f === 'mp4' || f === 'webm')));
        const hasSamples = api.options() && 'samples' in api.options();
        show(samples, hasSamples);
        show(diskMotion, hasDisk());
        const base = fileName.value.replace(/\.(mp4|webm|zip|mpeg)$/i, '');
        fileName.value = base + extension(f);
        const notes = [];
        if (viaHelper) {
          notes.push(`Frames go to the output folder through the helper (render_server.py)${f === 'h-png' ? '' : `, then ffmpeg makes ${base}${extension(f)} there`}.`);
          if (f === 'h-mpeg' && !MPEG_RATES.some((r) => Math.abs(r - Number(fps.value)) < 1e-3)) notes.push(`MPEG-2 allows only ${MPEG_RATES.join(', ')} fps; choose one of them.`);
          if (helper && !helper.ffmpeg && f !== 'h-png') notes.push('The helper did not find ffmpeg on its PATH.');
        } else {
          notes.push(fileHandle ? `Saved to ${fileHandle.name}.` : 'Saved through the browser\'s download: it asks for the place if it is set to (Firefox: Settings > Downloads > "Always ask you where to save files"; Safari: Settings > General > File download location > Ask for each download).');
          if (!helper) notes.push('To write into a folder you type here, and to make .mpeg files with ffmpeg, open the viewer from the helper: python3 render_server.py (in the viewer folder), then http://127.0.0.1:8765/.');
          if (f === 'zip') notes.push('The whole .zip is assembled in memory before it is saved.');
        }
        notes.push('The camera stays where it is now; every frame is traced until the style reports it final.');
        command.value = commandLine();
        destNote.textContent = notes.join('\n');
        destNote.classList.toggle('warn', notes.some((n) => /only|did not/.test(n)));
      };

      async function open() {
        if (!api.hasScene() || job) return;
        helper = await helperInfo();
        resolution.textContent = '';
        const view = api.view;
        const vw = even(view.clientWidth * (window.devicePixelRatio || 1)), vh = even(view.clientHeight * (window.devicePixelRatio || 1));
        for (const [w, h] of RESOLUTIONS) resolution.append(el('option', { value: `${w}x${h}`, textContent: `${w} × ${h}` }));
        resolution.append(el('option', { value: `${vw}x${vh}`, textContent: `this view (${vw} × ${vh})` }), el('option', { value: 'custom', textContent: 'custom' }));
        resolution.value = resolution.dataset.chosen || '1920x1080';
        if (!resolution.value) resolution.value = '1920x1080';
        format.textContent = '';
        for (const f of formats()) format.append(el('option', { value: f.value, textContent: f.label }));
        if (!format.dataset.chosen || ![...format.options].some((o) => o.value === format.dataset.chosen)) format.value = helper ? 'h-mpeg' : (typeof VideoEncoder === 'function' ? 'mp4' : 'zip');
        else format.value = format.dataset.chosen;
        const [a, b] = api.timeRange();
        t0.value = String(+a.toPrecision(8)); t1.value = String(+b.toPrecision(8));
        if (hasDisk()) diskMotion.checked = !!api.options().diskMotion;
        // The style's suggested slow-motion factor, until one is typed here.
        const win = slowWindow();
        if (win && !slow.dataset.chosen) slow.value = win.factor ?? 4;
        fileName.value = safeName(`${api.sceneName()}_${api.styleId()}`);
        if (helper && !folder.value) folder.value = helper.defaultFolder || '~/Documents/Visualization/renders';
        fileHandle = null;
        progress.classList.add('hidden');
        status.textContent = '';
        start.disabled = false;
        close.textContent = 'Close';
        refresh();
        dialog.showModal();
      }
      button.addEventListener('click', open);
      [resolution, width, height, fps, t0, t1, seconds, format, slow, samples, folder, fileName, keepFrames, diskMotion].forEach((i) => i.addEventListener('input', () => {
        if (i === resolution) resolution.dataset.chosen = resolution.value;
        if (i === slow) slow.dataset.chosen = '1';
        if (i === format) { format.dataset.chosen = format.value; fileHandle = null; }
        refresh();
      }));
      chooseFile.addEventListener('click', async () => {
        try {
          const mp4 = format.value === 'mp4';
          fileHandle = await window.showSaveFilePicker({ suggestedName: fileName.value, types: [{ description: mp4 ? 'MP4 video' : 'WebM video', accept: mp4 ? { 'video/mp4': ['.mp4'] } : { 'video/webm': ['.webm'] } }] });
          fileName.value = fileHandle.name;
          refresh();
        } catch (e) { /* cancelled */ }
      });
      close.addEventListener('click', () => {
        if (job) { job.cancelled = true; close.textContent = 'Cancelling…'; return; }
        dialog.close();
      });
      dialog.addEventListener('cancel', (e) => { if (job) { e.preventDefault(); job.cancelled = true; } });

      // Draw frame k until the style reports it final, then capture it in the same task (the drawing buffer is
      // still valid there).
      const renderFrameFinal = (cfg, capture) => new Promise((resolve, reject) => {
        let n = 0;
        const step = () => {
          if (job.cancelled) { resolve(false); return; }
          try {
            api.renderFrame(n === 0 ? 1 / cfg.fps : 0);
            n++;
            const style = api.style();
            const done = style && style.converged ? style.converged() : n >= 2;
            if (done) { capture(api.renderer.domElement); resolve(true); return; }
          } catch (e) { reject(e); return; }
          requestAnimationFrame(step);
        };
        requestAnimationFrame(step);
      });

      // The render itself: frames firstFrame .. lastFrame - 1 of cfg.frames (all by default), each drawn until the
      // style reports it final and handed to the sink; report(progress) after each. Returns the sink's message.
      async function runRender(cfg, report) {
        const options = api.options();
        const savedSamples = options && 'samples' in options ? options.samples : undefined;
        const savedMotion = hasDisk() ? options.diskMotion : undefined;
        let sink = null;
        const began = performance.now();
        try {
          sink = cfg.format.startsWith('h-') ? await helperSink(cfg, helper)
            : cfg.format === 'zip' ? zipSink(cfg) : await encoderSink(cfg, cfg.format);
          if (savedSamples !== undefined) options.samples = Math.max(1, Math.round(cfg.samples));
          if (savedMotion !== undefined) options.diskMotion = cfg.diskMotion ? 1 : 0;
          const times = frameTimes(cfg.t0, cfg.t1, cfg.frames, cfg.slow);
          window.GRRender.lastTimes = times;
          const first = Math.max(0, cfg.firstFrame || 0), last = Math.min(cfg.frames, cfg.lastFrame ?? cfg.frames);
          api.begin(cfg.width, cfg.height);
          for (let k = first; k < last; k++) {
            const t = times[k];
            api.setClock(t);
            const ok = await renderFrameFinal(cfg, (canvas) => sink.capture(canvas, k));
            if (!ok) break;
            await sink.flush();
            const elapsed = (performance.now() - began) / 1000, done = k - first + 1;
            report({ k, done, total: last - first, frames: cfg.frames, t, label: api.clockLabel(), describe: sink.describe, elapsed, left: elapsed / done * (last - first - done) });
          }
          if (job.cancelled) { sink.abort(); return null; }
          return await sink.finish();
        } catch (e) {
          if (sink) sink.abort();
          throw e;
        } finally {
          api.end();
          if (savedSamples !== undefined) options.samples = savedSamples;
          if (savedMotion !== undefined) options.diskMotion = savedMotion;
        }
      }

      start.addEventListener('click', async () => {
        const cfg = {
          width: even(Number(width.value)), height: even(Number(height.value)), fps: Number(fps.value),
          t0: Number(t0.value), t1: Number(t1.value), frames: frames(), mbps: Number(mbps.value), format: format.value,
          fileName: fileName.value || 'render', baseName: safeName(fileName.value.replace(/\.(mp4|webm|zip|mpeg)$/i, '')),
          folder: folder.value, keepFrames: keepFrames.checked, fileHandle, slow: slowWindow() ? Math.max(1, Number(slow.value) || 1) : 1,
          samples: Number(samples.value), diskMotion: diskMotion.checked,
        };
        if (!(cfg.fps > 0 && cfg.frames >= 1 && cfg.width >= 16 && cfg.height >= 16 && Number.isFinite(cfg.t0) && Number.isFinite(cfg.t1))) {
          status.textContent = 'Check the numbers.'; return;
        }
        job = { cancelled: false };
        start.disabled = true;
        close.textContent = 'Cancel';
        progress.classList.remove('hidden');
        progress.value = 0;
        const began = performance.now();
        try {
          status.textContent = 'Starting…';
          const message = await runRender(cfg, (p) => {
            progress.value = p.done / p.total;
            status.textContent = `${p.describe}\nframe ${p.k + 1} of ${p.frames}, ${p.label} = ${api.fmt(p.t)}\nelapsed ${clockText(p.elapsed)}, about ${clockText(p.left)} left`;
          });
          status.textContent = message === null ? 'Cancelled.' : `Done in ${clockText((performance.now() - began) / 1000)}: ${message}`;
        } catch (e) {
          status.textContent = 'Render failed: ' + (e && e.message ? e.message : e);
        } finally {
          job = null;
          start.disabled = false;
          close.textContent = 'Close';
        }
      });

      // Scripted renders (render_cli.py): the same render as the dialog's, from a settings object. Defaults: the
      // viewer's time range, 24 fps, 10 s, 1920 x 1080, 8 samples, the style's slow-motion factor, H.264 through
      // the helper. firstFrame / lastFrame render part of the frames (for splitting a render across machines).
      window.GRRender.render = async (o = {}) => {
        if (job) throw new Error('a render is already running');
        if (!api.hasScene()) throw new Error('no scene loaded');
        helper = await helperInfo();
        const [a, b] = api.timeRange(), win = slowWindow();
        const fpsV = o.fps ?? 24;
        const cfg = {
          width: even(o.width ?? 1920), height: even(o.height ?? 1080), fps: fpsV,
          t0: o.t0 ?? a, t1: o.t1 ?? b, frames: o.frames ?? Math.max(1, Math.round(fpsV * (o.seconds ?? 10))),
          mbps: o.mbps ?? 12, format: o.format ?? 'h-mp4', fileName: o.name ?? 'render', baseName: safeName(o.name ?? 'render'),
          folder: o.folder ?? '~/Documents/Visualization/renders', keepFrames: !!o.keepFrames, fileHandle: null,
          slow: win ? Math.max(1, o.slow ?? win.factor ?? 1) : 1, samples: o.samples ?? 8, diskMotion: !!o.diskMotion,
          firstFrame: o.firstFrame, lastFrame: o.lastFrame,
        };
        if (cfg.format.startsWith('h-') && !helper) throw new Error('this format needs the page to be served by render_server.py');
        job = { cancelled: false };
        try { return await runRender(cfg, o.onProgress || (() => {})); }
        finally { job = null; }
      };
      window.GRRender.cancel = () => { if (job) job.cancelled = true; };
      // For automated tests.
      window.GRRender.open = open;
      window.GRRender.dialog = dialog;
    },
  };
})();
