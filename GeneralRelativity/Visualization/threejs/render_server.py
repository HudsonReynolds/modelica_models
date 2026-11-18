#!/usr/bin/env python3
"""Local helper for the GeneralRelativity viewer's Render dialog.

Serves the viewer at http://127.0.0.1:8765/ and, when a render is started from that page, writes the frames
into a folder typed in the dialog and runs ffmpeg to make an MPEG-2 (.mpeg) or H.264 (.mp4) video there.

    python3 render_server.py                # serve the folder this script is in
    python3 render_server.py --port 9000 --root ~/Documents/Visualization

Safety: it listens on 127.0.0.1 only; every request that writes files must carry a random token that only
pages served by this helper can read (/api/info is same-origin), the Host header must be this server's
(no DNS rebinding), and output folders must be inside your home folder.
"""
import argparse
import http.server
import json
import os
import re
import secrets
import shutil
import subprocess
import sys
import urllib.parse

HOME = os.path.realpath(os.path.expanduser('~'))
TOKEN = secrets.token_urlsafe(24)
FFMPEG = shutil.which('ffmpeg')
# MPEG-2 frame rates as ffmpeg rationals.
MPEG_RATES = {'23.976': '24000/1001', '24': '24', '25': '25', '29.97': '30000/1001', '30': '30', '50': '50', '59.94': '60000/1001', '60': '60'}
job = None   # the render in progress: dict(folder, frames, name, format, fps, keep, count)


def rate_text(fps):
    text = ('%.3f' % fps).rstrip('0').rstrip('.')
    return text


def encode(frames, out, fps):
    """Make a video from frames/frame_%05d.png with ffmpeg: MPEG-2 for .mpeg/.mpg, else H.264. Returns
    (ok, error text)."""
    mpeg = out.lower().endswith(('.mpeg', '.mpg'))
    if mpeg and rate_text(fps) not in MPEG_RATES:
        return False, 'MPEG-2 allows only %s fps' % ', '.join(MPEG_RATES)
    rate = MPEG_RATES[rate_text(fps)] if mpeg else rate_text(fps)
    command = [FFMPEG, '-y', '-hide_banner', '-loglevel', 'error', '-framerate', rate,
               '-i', os.path.join(frames, 'frame_%05d.png')]
    if mpeg:
        # MPEG-2 video in an MPEG program stream; -q:v 2 is near the best quality the format allows.
        command += ['-c:v', 'mpeg2video', '-q:v', '2', '-pix_fmt', 'yuv420p', '-f', 'mpeg', out]
    else:
        command += ['-c:v', 'libx264', '-crf', '16', '-preset', 'slow', '-pix_fmt', 'yuv420p', '-movflags', '+faststart', out]
    result = subprocess.run(command, capture_output=True, text=True)
    return result.returncode == 0, result.stderr.strip()[-600:]


class Handler(http.server.SimpleHTTPRequestHandler):
    server_version = 'GRRenderHelper/1'

    def log_message(self, fmt, *args):
        if '/api/render/frame' not in str(args[0] if args else ''):
            sys.stderr.write('%s\n' % (fmt % args))

    def end_headers(self):
        # Revalidate every file, so an updated viewer is never run from a stale browser cache.
        if not self.path.startswith('/api/'):
            self.send_header('Cache-Control', 'no-cache')
        super().end_headers()

    # --- checks -----------------------------------------------------------------------------------------
    def host_ok(self):
        port = self.server.server_address[1]
        return self.headers.get('Host', '') in ('127.0.0.1:%d' % port, 'localhost:%d' % port)

    def origin_ok(self):
        origin = self.headers.get('Origin')
        port = self.server.server_address[1]
        return origin is None or origin in ('http://127.0.0.1:%d' % port, 'http://localhost:%d' % port)

    def reply(self, status, payload):
        body = json.dumps(payload).encode()
        self.send_response(status)
        self.send_header('Content-Type', 'application/json')
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        self.end_headers()
        self.wfile.write(body)

    def fail(self, status, message):
        self.reply(status, {'ok': False, 'error': message})

    # --- requests -----------------------------------------------------------------------------------------
    def do_GET(self):
        if not self.host_ok():
            return self.fail(403, 'unexpected Host header')
        if self.path.split('?')[0] == '/api/info':
            return self.reply(200, {'ok': True, 'token': TOKEN, 'ffmpeg': bool(FFMPEG),
                                    'defaultFolder': '~/Documents/Visualization/renders', 'root': os.path.realpath(self.directory)})
        return super().do_GET()

    def do_POST(self):
        global job
        if not (self.host_ok() and self.origin_ok()) or self.headers.get('X-Render-Token') != TOKEN:
            return self.fail(403, 'not allowed (open the viewer from this helper)')
        length = int(self.headers.get('Content-Length', 0))
        body = self.rfile.read(length) if length else b''
        url = urllib.parse.urlparse(self.path)
        try:
            if url.path == '/api/render/start':
                return self.api_start(json.loads(body or b'{}'))
            if url.path == '/api/render/frame':
                return self.api_frame(urllib.parse.parse_qs(url.query), body)
            if url.path == '/api/render/finish':
                return self.api_finish()
            if url.path == '/api/render/cancel':
                if job and not job['keep'] and os.path.isdir(job['frames']):
                    shutil.rmtree(job['frames'], ignore_errors=True)
                job = None
                return self.reply(200, {'ok': True})
        except Exception as e:   # report to the dialog rather than dropping the connection
            return self.fail(500, str(e))
        return self.fail(404, 'unknown request')

    # (Not named start/finish: StreamRequestHandler.finish() is called after every request.)
    def api_start(self, request):
        global job
        folder = os.path.realpath(os.path.expanduser(str(request.get('folder') or '~/Documents/Visualization/renders')))
        if folder != HOME and not folder.startswith(HOME + os.sep):
            return self.fail(400, 'the output folder must be inside your home folder (%s)' % HOME)
        name = re.sub(r'[^A-Za-z0-9._-]+', '_', str(request.get('name') or 'render')).strip('._') or 'render'
        fmt = request.get('format')
        if fmt not in ('h-mpeg', 'h-mp4', 'h-png'):
            return self.fail(400, 'unknown format')
        fps = float(request.get('fps') or 24)
        if fmt == 'h-mpeg' and rate_text(fps) not in MPEG_RATES:
            return self.fail(400, 'MPEG-2 allows only %s fps' % ', '.join(MPEG_RATES))
        if fmt != 'h-png' and not FFMPEG:
            return self.fail(400, 'ffmpeg was not found on the PATH; choose PNG frames or install ffmpeg')
        os.makedirs(folder, exist_ok=True)
        frames = os.path.join(folder, name if fmt == 'h-png' else name + '_frames')
        base, k = frames, 2
        while os.path.exists(frames) and os.listdir(frames):
            frames = '%s_%d' % (base, k)
            k += 1
        os.makedirs(frames, exist_ok=True)
        job = {'folder': folder, 'frames': frames, 'name': name, 'format': fmt, 'fps': fps,
               'keep': fmt == 'h-png' or bool(request.get('keepFrames')), 'count': 0}
        return self.reply(200, {'ok': True, 'folder': folder, 'frames': frames})

    def api_frame(self, query, body):
        if not job:
            return self.fail(409, 'no render started')
        i = int(query.get('i', ['-1'])[0])
        if not 0 <= i < 10 ** 6 or not body.startswith(b'\x89PNG'):
            return self.fail(400, 'bad frame')
        with open(os.path.join(job['frames'], 'frame_%05d.png' % i), 'wb') as f:
            f.write(body)
        job['count'] = max(job['count'], i + 1)
        return self.reply(200, {'ok': True})

    def api_finish(self):
        global job
        if not job:
            return self.fail(409, 'no render started')
        current, job = job, None
        if current['format'] == 'h-png':
            return self.reply(200, {'ok': True, 'message': '%d PNG frames in %s' % (current['count'], current['frames'])})
        mpeg = current['format'] == 'h-mpeg'
        out = os.path.join(current['folder'], current['name'] + ('.mpeg' if mpeg else '.mp4'))
        ok, error = encode(current['frames'], out, current['fps'])
        if not ok:
            return self.fail(500, 'ffmpeg failed (frames kept in %s): %s' % (current['frames'], error))
        if not current['keep']:
            shutil.rmtree(current['frames'], ignore_errors=True)
        size = os.path.getsize(out) / 1e6
        kept = '; frames kept in %s' % current['frames'] if current['keep'] else ''
        return self.reply(200, {'ok': True, 'message': 'wrote %s (%.1f MB, %d frames)%s' % (out, size, current['count'], kept), 'output': out})


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--port', type=int, default=8765)
    parser.add_argument('--root', default=os.path.dirname(os.path.abspath(__file__)), help='folder to serve (the viewer)')
    args = parser.parse_args()
    root = os.path.realpath(os.path.expanduser(args.root))
    handler = lambda *a, **k: Handler(*a, directory=root, **k)
    server = http.server.ThreadingHTTPServer(('127.0.0.1', args.port), handler)
    print('Viewer: http://127.0.0.1:%d/  (serving %s)' % (args.port, root))
    print('ffmpeg: %s' % (FFMPEG or 'not found: only PNG frames can be written to folders'))
    print('Press Ctrl+C to stop.')
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == '__main__':
    main()
