// The sound of a compact-binary merger: the gravitational-wave strain recorded by Examples.BinaryMerger
// (the h+ and hx scalars of the tracks, face-on, R h / M) shown as a "sound profile" (waveform and frequency
// against time, with a playhead on the clock) and played as audio.
//
// The physical time scale comes from the scene's gwSource record (seconds per unit M for the binary's total
// mass). For the default 36 + 29 solar-mass black holes the true sound is the GW150914-like chirp: about
// 17 Hz rising to a ~270 Hz ringdown in half a second. Playback options:
//   - real time: the strain itself, at its true frequencies and duration;
//   - slowed down (x4, x16) with the pitch kept, by resynthesis from the instantaneous amplitude and phase;
//   - pitch shift (x2 ... x8) for laptop speakers that cannot reproduce 20-60 Hz;
//   - a live tone that follows the animation clock (frequency and loudness of the wave at that moment);
//   - "Save .wav" writes the chosen version as a 16-bit WAV file.
// Classic script; the core calls GRAudio.install(api).
(function () {
  const RATE = 44100;
  const css = `
    #gwSound .row { display: flex; gap: 6px; align-items: center; flex-wrap: wrap; margin: 4px 0; font-size: 12px; }
    #gwSound select, #gwSound button { padding: 2px 6px; font-size: 12px; }
    #gwSound .note { color: #8a93a6; font-size: 11px; white-space: pre-line; }
  `;

  // Build the waveform from a track with h+ and hx scalars: amplitude, unwrapped GW phase and frequency.
  function analyse(data) {
    const spec = data.trackSpecs.find((t) => t.scalarNames.includes('h+') && t.scalarNames.includes('hx'));
    if (!spec) return null;
    const ip = spec.scalarNames.indexOf('h+'), ix = spec.scalarNames.indexOf('hx'), ifr = spec.scalarNames.indexOf('f_GW Hz');
    const c = data.header.c || 1;
    const n = spec.samples.length;
    const t = new Float64Array(n), hp = new Float64Array(n), amp = new Float64Array(n), phase = new Float64Array(n), freq = new Float64Array(n);
    let last = 0, offset = 0, peak = 0, iPeak = 0;
    spec.samples.forEach((s, k) => {
      t[k] = s.x[0] / c;
      hp[k] = s.scalars[ip];
      const re = s.scalars[ip], im = -s.scalars[ix];     // H = h+ - i hx ~ A e^{-i Phi}
      amp[k] = Math.hypot(re, im);
      let ph = -Math.atan2(im, re);
      if (k > 0) { while (ph + offset - last > Math.PI) offset -= 2 * Math.PI; while (ph + offset - last < -Math.PI) offset += 2 * Math.PI; }
      phase[k] = ph + offset;
      last = phase[k];
      if (amp[k] > peak) { peak = amp[k]; iPeak = k; }
    });
    for (let k = 0; k < n; k++) {
      const a = Math.max(0, k - 1), b = Math.min(n - 1, k + 1);
      freq[k] = (phase[b] - phase[a]) / Math.max(t[b] - t[a], 1e-30) / (2 * Math.PI);   // cycles per unit time (M)
    }
    // Where the ringdown has decayed to numerical noise, hold the last reliable frequency.
    for (let k = iPeak; k < n; k++) if (amp[k] < 1e-5 * peak) freq[k] = freq[k - 1];
    // Seconds per M: from the scene, else from the f_GW Hz scalar against the phase-derived frequency.
    const src = (data.records.gwSource || [])[0];
    let secondsPerM = src && src.secondsPerM;
    if (!secondsPerM && ifr >= 0) {
      const ratios = [];
      for (let k = 1; k < iPeak; k += Math.max(1, Math.floor(iPeak / 50))) { const hz = spec.samples[k].scalars[ifr]; if (hz > 0 && freq[k] > 0) ratios.push(freq[k] / hz); }
      ratios.sort((x, y) => x - y);
      secondsPerM = ratios.length ? ratios[ratios.length >> 1] : null;
    }
    const metric = data.metric && data.metric.name === 'binary' ? data.metric : null;
    const tMerge = metric && metric.p[4] < 1e29 ? metric.p[4] / c : t[iPeak];
    return { spec, n, t, hp, amp, phase, freq, peak, iPeak, secondsPerM, tMerge, massSolar: src && src.totalMassSolar, ringdownHz: src && src.ringdownHz };
  }

  // Interpolation helpers on the sample times.
  function locate(w, time) {
    const t = w.t, n = w.n;
    if (time <= t[0]) return [0, 0];
    if (time >= t[n - 1]) return [n - 2, 1];
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (t[mid] <= time) lo = mid; else hi = mid; }
    return [lo, (time - t[lo]) / Math.max(t[lo + 1] - t[lo], 1e-30)];
  }
  const lerp = (arr, [i, f]) => arr[i] + f * (arr[i + 1] - arr[i]);

  // Audio samples: the strain resynthesized from amplitude and frequency, time stretched by `stretch` (pitch
  // kept) and the pitch multiplied by `pitch`; normalized to the peak, with short fades.
  function synthesize(w, stretch, pitch) {
    const T = (w.t[w.n - 1] - w.t[0]) * w.secondsPerM * stretch;
    const N = Math.max(1, Math.round(T * RATE));
    const out = new Float32Array(N);
    let ph = 0;
    for (let s = 0; s < N; s++) {
      const time = w.t[0] + (s / RATE) / (w.secondsPerM * stretch);
      const at = locate(w, time);
      const hz = Math.max(0, lerp(w.freq, at)) / w.secondsPerM * pitch;
      ph += 2 * Math.PI * Math.min(hz, 0.45 * RATE) / RATE;
      const fade = Math.min(1, s / (0.02 * RATE), (N - 1 - s) / (0.02 * RATE));
      out[s] = 0.9 * (lerp(w.amp, at) / w.peak) * Math.sin(ph) * fade;
    }
    return out;
  }
  function wav(samples) {
    const buf = new ArrayBuffer(44 + 2 * samples.length), v = new DataView(buf);
    const str = (o, s) => [...s].forEach((ch, i) => v.setUint8(o + i, ch.charCodeAt(0)));
    str(0, 'RIFF'); v.setUint32(4, 36 + 2 * samples.length, true); str(8, 'WAVE'); str(12, 'fmt ');
    v.setUint32(16, 16, true); v.setUint16(20, 1, true); v.setUint16(22, 1, true); v.setUint32(24, RATE, true);
    v.setUint32(28, 2 * RATE, true); v.setUint16(32, 2, true); v.setUint16(34, 16, true); str(36, 'data'); v.setUint32(40, 2 * samples.length, true);
    for (let i = 0; i < samples.length; i++) v.setInt16(44 + 2 * i, Math.max(-1, Math.min(1, samples[i])) * 32767, true);
    return new Blob([buf], { type: 'audio/wav' });
  }

  window.GRAudio = {
    install(api) {
      document.head.append(Object.assign(document.createElement('style'), { textContent: css }));
      const section = document.createElement('div');
      section.id = 'gwSound';
      section.style.display = 'none';
      section.innerHTML = `<h2>Gravitational-wave sound</h2>
        <canvas width="320" height="170"></canvas>
        <div class="row"><button type="button" data-act="play">▶ Play chirp</button>
          <select data-opt="stretch" title="Slow the sound down; the pitch is kept">
            <option value="1">real time</option><option value="4">4× slower</option><option value="16">16× slower</option></select>
          <select data-opt="pitch" title="Raise the pitch (laptop speakers cannot play the lowest frequencies)">
            <option value="1">true pitch</option><option value="2">pitch ×2</option><option value="4">pitch ×4</option><option value="8">pitch ×8</option></select>
          <button type="button" data-act="wav" title="Save the chosen version as a WAV file">Save .wav</button></div>
        <div class="row"><label><input type="checkbox" data-opt="live"> live tone with the animation</label></div>
        <div class="note"></div>`;
      const legendTitle = [...api.aside.querySelectorAll('h2')].find((h) => h.textContent === 'Bodies');
      api.aside.insertBefore(section, legendTitle || null);
      const canvas = section.querySelector('canvas'), note = section.querySelector('.note');
      const q = (sel) => section.querySelector(sel);
      let wave = null, backdrop = null, ctx = null, osc = null, gain = null, source = null, audioStart = 0, audioStretch = 1;

      const audio = () => {
        if (!ctx) ctx = new (window.AudioContext || window.webkitAudioContext)();
        if (ctx.state === 'suspended') ctx.resume();
        return ctx;
      };
      const seconds = (time) => (time - wave.tMerge) * wave.secondsPerM;

      // The sound profile: h+ (normalized) and the frequency (log scale) against time to the merger.
      function drawBackdrop() {
        const W = canvas.width, H = canvas.height, g = canvas.getContext('2d');
        g.clearRect(0, 0, W, H);
        g.fillStyle = '#0f131b'; g.fillRect(0, 0, W, H);
        const t0 = wave.t[0], t1 = wave.t[wave.n - 1];
        const X = (time) => 30 + (W - 36) * (time - t0) / (t1 - t0);
        const top = 8, mid = 58, half = 44;
        // Waveform.
        g.strokeStyle = '#5fb3ff'; g.lineWidth = 1; g.beginPath();
        const step = Math.max(1, Math.floor(wave.n / 2000));
        for (let k = 0; k < wave.n; k += step) {
          const x = X(wave.t[k]), y = mid - half * wave.hp[k] / wave.peak;
          if (k === 0) g.moveTo(x, y); else g.lineTo(x, y);
        }
        g.stroke();
        // Frequency on a log axis below.
        const fTop = 112, fH = 40;
        const hz = (k) => wave.freq[k] / wave.secondsPerM;
        let fMin = Infinity, fMax = 0;
        for (let k = 0; k < wave.n; k++) { const f = hz(k); if (f > 0) { fMin = Math.min(fMin, f); fMax = Math.max(fMax, f); } }
        fMin = Math.max(fMin, 1); fMax = Math.max(fMax, fMin * 2);
        const Y = (f) => fTop + fH - fH * (Math.log(Math.max(f, fMin)) - Math.log(fMin)) / (Math.log(fMax) - Math.log(fMin));
        g.strokeStyle = '#ffb454'; g.beginPath();
        for (let k = 0; k < wave.n; k += step) { const x = X(wave.t[k]), y = Y(hz(k)); if (k === 0) g.moveTo(x, y); else g.lineTo(x, y); }
        g.stroke();
        // Merger marker and labels.
        g.strokeStyle = 'rgba(255,255,255,0.35)'; g.setLineDash([3, 3]); g.beginPath();
        g.moveTo(X(wave.tMerge), top); g.lineTo(X(wave.tMerge), fTop + fH); g.stroke(); g.setLineDash([]);
        g.fillStyle = '#8a93a6'; g.font = '10px system-ui';
        g.fillText('h+', 4, mid + 3);
        g.fillText(`${Math.round(fMax)} Hz`, 2, fTop + 8);
        g.fillText(`${Math.round(fMin)} Hz`, 2, fTop + fH);
        g.fillText(`${(seconds(t0)).toFixed(seconds(t0) < -2 ? 1 : 2)} s`, 30, H - 4);
        g.fillText('merger', Math.min(X(wave.tMerge) - 34, W - 120), top + 8);
        const end = `+${(1000 * seconds(t1)).toFixed(0)} ms`;
        g.fillText(end, W - 4 - g.measureText(end).width, H - 4);
        backdrop = g.getImageData(0, 0, W, H);
        return X;
      }
      let X = null;
      function drawPlayheads(clock) {
        if (!wave || !backdrop) return;
        const g = canvas.getContext('2d');
        g.putImageData(backdrop, 0, 0);
        const line = (time, colour) => {
          if (time < wave.t[0] || time > wave.t[wave.n - 1]) return;
          g.strokeStyle = colour; g.lineWidth = 1.5; g.beginPath(); g.moveTo(X(time), 6); g.lineTo(X(time), canvas.height - 14); g.stroke();
        };
        line(clock, '#e6e9ef');
        if (source && ctx) line(wave.t[0] + (ctx.currentTime - audioStart) / (wave.secondsPerM * audioStretch), '#7ee0a1');
      }

      api.onScene((data) => {
        stopAll();
        wave = analyse(data);
        section.style.display = wave && wave.secondsPerM ? '' : 'none';
        if (!wave || !wave.secondsPerM) return;
        X = drawBackdrop();
        const dur = (wave.t[wave.n - 1] - wave.t[0]) * wave.secondsPerM;
        const f0 = wave.freq[0] / wave.secondsPerM;
        const before = locate(wave, wave.tMerge - 1e-6 * Math.abs(wave.tMerge));
        const fMerge = lerp(wave.freq, before) / wave.secondsPerM;
        const ring = wave.ringdownHz ? `, then a ${wave.ringdownHz.toFixed(0)} Hz ringdown` : '';
        note.textContent = `${wave.massSolar ? wave.massSolar + ' solar masses: ' : ''}${f0.toFixed(f0 < 100 ? 1 : 0)} Hz rising to ${fMerge.toFixed(0)} Hz at the merger${ring}; ${dur < 1 ? (dur * 1000).toFixed(0) + ' ms' : dur.toFixed(2) + ' s'} in real time. `
          + 'Face-on strain h+ from the simulation; the sound is resynthesized from its amplitude and phase. Below about 60 Hz most laptop speakers are silent: raise the pitch or use headphones.';
      });

      function stopAll() {
        if (source) { try { source.stop(); } catch (e) { /* stopped */ } source = null; }
        q('[data-act=play]').textContent = '▶ Play chirp';
      }
      q('[data-act=play]').addEventListener('click', () => {
        if (!wave) return;
        if (source) { stopAll(); return; }
        const a = audio();
        audioStretch = Number(q('[data-opt=stretch]').value);
        const samples = synthesize(wave, audioStretch, Number(q('[data-opt=pitch]').value));
        const buffer = a.createBuffer(1, samples.length, RATE);
        buffer.copyToChannel(samples, 0);
        source = a.createBufferSource();
        source.buffer = buffer;
        source.connect(a.destination);
        source.onended = () => { source = null; q('[data-act=play]').textContent = '▶ Play chirp'; };
        audioStart = a.currentTime;
        source.start();
        q('[data-act=play]').textContent = '■ Stop';
      });
      q('[data-act=wav]').addEventListener('click', () => {
        if (!wave) return;
        const stretch = Number(q('[data-opt=stretch]').value), pitch = Number(q('[data-opt=pitch]').value);
        const name = `${(api.data().header.name || 'merger').replace(/[^A-Za-z0-9._-]+/g, '_')}_chirp${stretch > 1 ? '_slow' + stretch : ''}${pitch > 1 ? '_pitch' + pitch : ''}.wav`;
        const a = Object.assign(document.createElement('a'), { href: URL.createObjectURL(wav(synthesize(wave, stretch, pitch))), download: name });
        document.body.append(a); a.click(); a.remove();
        setTimeout(() => URL.revokeObjectURL(a.href), 60000);
      });
      q('[data-opt=live]').addEventListener('change', (e) => {
        if (e.target.checked) {
          const a = audio();
          osc = a.createOscillator(); gain = a.createGain();
          gain.gain.value = 0; osc.connect(gain); gain.connect(a.destination); osc.start();
        } else if (osc) { osc.stop(); osc = null; gain = null; }
      });

      // Every frame: the playhead, and the live tone at the wave's frequency and loudness at the clock.
      api.onFrame((clock, playing) => {
        if (!wave || section.style.display === 'none') return;
        drawPlayheads(clock);
        if (osc && ctx) {
          const inside = clock >= wave.t[0] && clock <= wave.t[wave.n - 1];
          const at = locate(wave, clock);
          const hz = Math.min(Math.max(lerp(wave.freq, at) / wave.secondsPerM * Number(q('[data-opt=pitch]').value), 1), 20000);
          const loud = playing && inside ? 0.35 * lerp(wave.amp, at) / wave.peak : 0;
          osc.frequency.setTargetAtTime(hz, ctx.currentTime, 0.015);
          gain.gain.setTargetAtTime(loud, ctx.currentTime, 0.03);
        }
      });
      window.GRAudio.analyse = analyse;      // for tests
      window.GRAudio.synthesize = synthesize;
      window.GRAudio.wav = wav;
    },
  };
})();
