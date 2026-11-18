#!/usr/bin/env python3
"""Render a viewer animation from a terminal, with the same renderer as the viewer's Render dialog.

It serves this folder on 127.0.0.1, opens the viewer in a headless browser (Firefox, else Chrome or
Chromium), loads the scene, applies the style, options and camera, traces every frame on the GPU exactly as
the Render dialog does, and writes PNG frames and a video (ffmpeg) through render_server.py. The browser does
not need to be open, and the render does not pause when you switch windows.

    python3 render_cli.py binary_neutron_star_merger.ndjson -o ~/Documents/Visualization/renders/merger.mp4
    python3 render_cli.py kerr_scene.ndjson --style raytrace --seconds 20 --size 3840x2160 -o kerr.mp4

The Render dialog's "command line" row gives the command for what is set up in the viewer (camera, options,
time range, frames). Output: .mp4 (H.264) or .mpeg (MPEG-2) by its extension, or a folder (a path ending in /)
for PNG frames.

Splitting a render across machines (frames are independent and reproducible): render parts of the frames as
PNG into one folder, then join them:

    python3 render_cli.py scene.ndjson <same settings> --frames 0:120 -o frames/      (machine 1)
    python3 render_cli.py scene.ndjson <same settings> --frames 120:240 -o frames/    (machine 2)
    python3 render_cli.py --encode frames/ --fps 24 -o merger.mp4

The GPU does the work, so a machine with a faster GPU renders faster. The "WebGL:" line names the renderer the
browser uses; a software one (SwiftShader, llvmpipe) is very slow.
"""
import argparse
import glob
import http.server
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import urllib.parse

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import render_server as rs  # noqa: E402  (the helper: token checks, frame upload, ffmpeg)

# Runs in the page: waits for the viewer, loads the scene, applies the settings and calls GRRender.render().
DRIVER = r"""
(async () => {
  const info = await (await fetch('/api/info', { cache: 'no-store' })).json();
  const post = (path, body) => fetch(path, { method: 'POST', headers: { 'Content-Type': 'application/json', 'X-Render-Token': info.token }, body: JSON.stringify(body) });
  try {
    for (let k = 0; k < 600 && !(window.GRViewer && window.GRRender && window.GRRender.render); k++) await new Promise((r) => setTimeout(r, 50));
    if (!(window.GRViewer && window.GRRender && window.GRRender.render)) throw new Error('the viewer did not start (it needs WebGL 2)');
    const cfg = await (await fetch('/cli/config', { cache: 'no-store' })).json();
    const V = window.GRViewer;
    const gl = V.renderer.getContext(), dbg = gl.getExtension('WEBGL_debug_renderer_info');
    await post('/cli/log', { message: 'WebGL: ' + gl.getParameter(dbg ? dbg.UNMASKED_RENDERER_WEBGL : gl.RENDERER) });
    V.loadScene(await (await fetch('/cli/scene', { cache: 'no-store' })).text());
    if (cfg.style) {
      const known = [...document.querySelectorAll('#styleSelect option')].filter((o) => !o.disabled).map((o) => o.value);
      if (!known.includes(cfg.style)) throw new Error(`style "${cfg.style}" is not available for this scene (available: ${known.join(', ')})`);
      V.activateStyle(cfg.style);
    }
    for (const [k, v] of Object.entries(cfg.options || {})) V.setOption(k, v);
    V.setPlaying(false);
    for (const name of cfg.hide || []) {
      const box = [...document.querySelectorAll('#legend .track')].find((e) => e.querySelector('.track-title') && e.querySelector('.track-title').textContent.trim() === name);
      const input = box && box.querySelector('input[type=checkbox]');
      if (!input) throw new Error(`no body named "${name}"`);
      if (input.checked) { input.checked = false; input.dispatchEvent(new Event('change')); }
    }
    if (cfg.noLegend) document.getElementById('showLegend').checked = false;
    if (cfg.camera) {
      V.camera.position.set(...cfg.camera);
      if (cfg.target) V.controls.target.set(...cfg.target);
      V.controls.update();
    }
    const message = await window.GRRender.render({
      width: cfg.width, height: cfg.height, fps: cfg.fps, seconds: cfg.seconds, t0: cfg.from ?? undefined, t1: cfg.to ?? undefined,
      samples: cfg.samples, slow: cfg.slow ?? undefined, diskMotion: cfg.diskRotates, format: cfg.format,
      folder: cfg.folder, name: cfg.name, keepFrames: cfg.keepFrames, firstFrame: cfg.firstFrame ?? undefined, lastFrame: cfg.lastFrame ?? undefined,
      onProgress: (p) => post('/cli/progress', p),
    });
    await post('/cli/done', { ok: true, message });
  } catch (e) {
    await post('/cli/done', { ok: false, error: String(e && e.message ? e.message : e) });
  }
})();
"""


def clock_text(s):
    m = int(s // 60)
    return '%d min %02d s' % (m, round(s - 60 * m)) if m else '%d s' % round(s)


class State:
    config = {}
    scene = b''
    job = None
    done = threading.Event()
    result = None
    started = False
    written = 0


class CLIHandler(rs.Handler):
    def log_message(self, fmt, *args):
        pass

    def token_ok(self):
        return self.host_ok() and self.origin_ok() and self.headers.get('X-Render-Token') == rs.TOKEN

    def do_GET(self):
        if not self.host_ok():
            return self.fail(403, 'unexpected Host header')
        path = self.path.split('?')[0]
        if path == '/cli.html':
            with open(os.path.join(HERE, 'index.html'), encoding='utf-8') as f:
                html = f.read()
            html = html.replace('</body>', '<script src="/cli/driver.js"></script>\n</body>', 1)
            return self.send_bytes(html.encode(), 'text/html; charset=utf-8')
        if path == '/cli/driver.js':
            return self.send_bytes(DRIVER.encode(), 'text/javascript; charset=utf-8')
        if path == '/cli/config':
            return self.send_bytes(json.dumps(State.config).encode(), 'application/json')
        if path == '/cli/scene':
            return self.send_bytes(State.scene, 'application/x-ndjson')
        return super().do_GET()

    def send_bytes(self, body, kind):
        self.send_response(200)
        self.send_header('Content-Type', kind)
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        self.end_headers()
        self.wfile.write(body)

    def do_POST(self):
        path = urllib.parse.urlparse(self.path).path
        if path.startswith('/cli/'):
            if not self.token_ok():
                return self.fail(403, 'not allowed')
            length = int(self.headers.get('Content-Length', 0))
            data = json.loads(self.rfile.read(length) or b'{}')
            if path == '/cli/log':
                State.started = True
                print(data.get('message', ''), flush=True)
            elif path == '/cli/progress':
                State.started = True
                line = 'frame %d of %d (%s = %.6g)  elapsed %s, about %s left' % (
                    data['k'] + 1, data['frames'], data.get('label', 't'), data['t'], clock_text(data['elapsed']), clock_text(data['left']))
                if sys.stdout.isatty():
                    print('\r' + line.ljust(90), end='', flush=True)
                else:
                    print(line, flush=True)
            elif path == '/cli/done':
                State.result = data
                State.done.set()
            return self.reply(200, {'ok': True})
        return super().do_POST()

    def api_frame(self, query, body):
        State.written += 1
        return super().api_frame(query, body)

    # The job is set by the command line, not by the page: exact folders, and parts of a render may share one.
    def api_start(self, request):
        rs.job = dict(State.job)
        os.makedirs(rs.job['frames'], exist_ok=True)
        return self.reply(200, {'ok': True, 'folder': rs.job['folder'], 'frames': rs.job['frames']})


def find_browser(choice):
    firefox = ['/Applications/Firefox.app/Contents/MacOS/firefox', shutil.which('firefox'),
               r'C:\Program Files\Mozilla Firefox\firefox.exe']
    chrome = ['/Applications/Google Chrome.app/Contents/MacOS/Google Chrome', '/Applications/Chromium.app/Contents/MacOS/Chromium',
              shutil.which('google-chrome'), shutil.which('google-chrome-stable'), shutil.which('chromium'), shutil.which('chromium-browser'),
              r'C:\Program Files\Google\Chrome\Application\chrome.exe']
    if choice not in ('auto', 'firefox', 'chrome'):
        return ('chrome' if 'chrom' in os.path.basename(choice).lower() else 'firefox'), choice
    for kind, paths in (('firefox', firefox), ('chrome', chrome)):
        if choice in ('auto', kind):
            for p in paths:
                if p and os.path.exists(p):
                    return kind, p
    return None, None


def launch(kind, binary, url, profile, extra):
    if kind == 'firefox':
        with open(os.path.join(profile, 'user.js'), 'w') as f:
            for k, v in (('browser.shell.checkDefaultBrowser', 'false'), ('datareporting.policy.dataSubmissionEnabled', 'false'),
                         ('toolkit.telemetry.reportingpolicy.firstRun', 'false'), ('browser.sessionstore.resume_from_crash', 'false'),
                         ('webgl.force-enabled', 'true')):
                f.write('user_pref("%s", %s);\n' % (k, v))
        args = [binary, '--headless', '--no-remote', '--profile', profile, '--window-size=1280,800']
    else:
        angle = {'Darwin': ['--use-angle=metal'], 'Linux': ['--use-angle=vulkan', '--enable-features=Vulkan'], 'Windows': ['--use-angle=d3d11']}
        args = [binary, '--headless=new', '--user-data-dir=' + profile, '--no-first-run', '--no-default-browser-check',
                '--window-size=1280,800', '--ignore-gpu-blocklist', '--enable-gpu'] + angle.get(platform.system(), [])
    return subprocess.Popen(args + list(extra) + [url], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def parse_triple(text, what):
    try:
        v = [float(x) for x in text.split(',')]
        assert len(v) == 3
        return v
    except Exception:
        sys.exit('%s must be x,y,z' % what)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('scene', nargs='?', help='scene file (.ndjson)')
    ap.add_argument('-o', '--output', help='.mp4 (H.264), .mpeg (MPEG-2), or a folder ending in / for PNG frames; '
                    'default: renders/<scene>_<style>.mp4 next to the scene')
    ap.add_argument('--style', help="the viewer style (default: the scene's own)")
    ap.add_argument('--option', action='append', default=[], metavar='KEY=VALUE', help='a style option (repeatable)')
    ap.add_argument('--hide', action='append', default=[], metavar='NAME', help='hide a body (repeatable)')
    ap.add_argument('--no-legend', action='store_true', help='no legend in the picture')
    ap.add_argument('--camera', help='camera position x,y,z (viewer coordinates); default: the style\'s camera')
    ap.add_argument('--target', help='point the camera looks at, x,y,z (default 0,0,0 with --camera)')
    ap.add_argument('--from', dest='t_from', type=float, help="first clock value (default: the viewer's range)")
    ap.add_argument('--to', dest='t_to', type=float, help="last clock value (default: the viewer's range)")
    ap.add_argument('--fps', type=float, default=24.0)
    ap.add_argument('--seconds', type=float, default=10.0, help='video length')
    ap.add_argument('--size', default='1920x1080', help='WIDTHxHEIGHT')
    ap.add_argument('--samples', type=int, default=8, help='ray tracer samples per frame (8 is indistinguishable from 48 in a video)')
    ap.add_argument('--slow', type=float, help="slow motion around a merger (default: the style's, 4 for neutron stars)")
    ap.add_argument('--disk-rotates', action='store_true', help='the accretion disk pattern turns with the clock')
    ap.add_argument('--frames', help='render only frames A:B (0-based, B excluded) of the whole sequence, as PNG frames')
    ap.add_argument('--keep-frames', action='store_true', help='keep the PNG frames next to the video')
    ap.add_argument('--encode', metavar='FOLDER', help='only join the PNG frames in FOLDER into the video -o')
    ap.add_argument('--browser', default='auto', help='auto, firefox, chrome, or the path of a browser')
    ap.add_argument('--browser-arg', action='append', default=[], help='extra browser argument (repeatable)')
    ap.add_argument('--port', type=int, default=0, help='port to serve on (default: any free one)')
    args = ap.parse_args()

    if args.encode:
        if not args.output:
            sys.exit('--encode needs -o (the video)')
        if not rs.FFMPEG:
            sys.exit('ffmpeg was not found on the PATH')
        frames = os.path.realpath(os.path.expanduser(args.encode))
        n = len(glob.glob(os.path.join(frames, 'frame_*.png')))
        out = os.path.realpath(os.path.expanduser(args.output))
        ok, error = rs.encode(frames, out, args.fps)
        if not ok:
            sys.exit('ffmpeg failed: ' + error)
        print('wrote %s (%.1f MB, %d frames)' % (out, os.path.getsize(out) / 1e6, n))
        return

    if not args.scene:
        ap.error('a scene file is needed')
    scene = os.path.realpath(os.path.expanduser(args.scene))
    if not os.path.isfile(scene):
        sys.exit('no scene file %s' % args.scene)
    with open(scene, 'rb') as f:
        State.scene = f.read()
    try:
        width, height = (int(v) for v in args.size.lower().split('x'))
    except ValueError:
        sys.exit('--size must be WIDTHxHEIGHT')
    first = last = None
    if args.frames:
        try:
            a, b = args.frames.split(':')
            first, last = int(a), int(b)
            assert 0 <= first < last
        except Exception:
            sys.exit('--frames must be A:B with 0 <= A < B')

    # Output: a video file, or a folder of PNG frames.
    stem = os.path.splitext(os.path.basename(scene))[0]
    output = args.output or os.path.join(os.path.dirname(scene), 'renders', '%s_%s.mp4' % (stem, args.style or 'render'))
    output = os.path.expanduser(output)
    ext = os.path.splitext(output)[1].lower()
    png = output.endswith(('/', os.sep)) or os.path.isdir(output) or ext == ''
    if args.frames and not png:
        sys.exit('--frames writes PNG frames: give a folder ending in / as -o, then join the parts with --encode')
    if not png and ext not in ('.mp4', '.mpeg', '.mpg'):
        sys.exit('-o must end in .mp4, .mpeg or / (a folder for PNG frames)')
    if not png and not rs.FFMPEG:
        sys.exit('ffmpeg was not found on the PATH: install it, or give a folder ending in / for PNG frames')
    if not png and ext in ('.mpeg', '.mpg') and rs.rate_text(args.fps) not in rs.MPEG_RATES:
        sys.exit('MPEG-2 allows only %s fps' % ', '.join(rs.MPEG_RATES))
    output = os.path.realpath(output)
    if png:
        folder, name, frames, fmt = output, os.path.basename(output), output, 'h-png'
    else:
        folder, name = os.path.dirname(output), os.path.splitext(os.path.basename(output))[0]
        frames, fmt = os.path.join(folder, name + '_frames'), 'h-mpeg' if ext in ('.mpeg', '.mpg') else 'h-mp4'
        # Stale frames of an earlier, longer render would end up in the video.
        for old in glob.glob(os.path.join(frames, 'frame_*.png')):
            os.remove(old)
    os.makedirs(frames, exist_ok=True)
    State.job = {'folder': folder, 'frames': frames, 'name': name, 'format': fmt, 'fps': args.fps,
                 'keep': png or args.keep_frames, 'count': 0}

    options = {}
    for item in args.option:
        if '=' not in item:
            sys.exit('--option must be KEY=VALUE')
        k, v = item.split('=', 1)
        try:
            options[k] = float(v)
        except ValueError:
            options[k] = v
    camera = parse_triple(args.camera, '--camera') if args.camera else None
    target = parse_triple(args.target, '--target') if args.target else ([0.0, 0.0, 0.0] if camera else None)
    State.config = {
        'style': args.style, 'options': options, 'hide': args.hide, 'noLegend': args.no_legend,
        'camera': camera, 'target': target, 'from': args.t_from, 'to': args.t_to, 'fps': args.fps, 'seconds': args.seconds,
        'width': width, 'height': height, 'samples': args.samples, 'slow': args.slow, 'diskRotates': args.disk_rotates,
        'format': fmt, 'folder': folder, 'name': name, 'keepFrames': args.keep_frames, 'firstFrame': first, 'lastFrame': last,
    }

    kind, binary = find_browser(args.browser)
    if not binary:
        sys.exit('no browser found: install Firefox or Chrome, or give --browser PATH')
    server = http.server.ThreadingHTTPServer(('127.0.0.1', args.port), lambda *a, **k: CLIHandler(*a, directory=HERE, **k))
    port = server.server_address[1]
    threading.Thread(target=server.serve_forever, daemon=True).start()
    profile = tempfile.mkdtemp(prefix='gr_render_')
    print('%s -> %s (%s, headless %s)' % (os.path.basename(scene), output + ('' if not png else os.sep), fmt[2:], kind), flush=True)
    proc = launch(kind, binary, 'http://127.0.0.1:%d/cli.html' % port, profile, args.browser_arg)
    began = time.time()
    try:
        while not State.done.wait(0.5):
            if proc.poll() is not None:
                sys.exit('\nthe browser exited before the render finished')
            if not State.started and time.time() - began > 120:
                sys.exit('\nthe page did not start within 2 minutes')
        result = State.result
        if sys.stdout.isatty():
            print()
        if not result.get('ok'):
            sys.exit('render failed: %s' % result.get('error'))
        message = result.get('message')
        if args.frames:
            message = '%d PNG frames (%d to %d of the sequence) in %s' % (State.written, first, first + State.written - 1, frames)
        print('done in %s: %s' % (clock_text(time.time() - began), message))
    except KeyboardInterrupt:
        print('\ncancelled')
        sys.exit(130)
    finally:
        proc.terminate()
        try:
            proc.wait(10)
        except subprocess.TimeoutExpired:
            proc.kill()
        server.shutdown()
        shutil.rmtree(profile, ignore_errors=True)


if __name__ == '__main__':
    main()
