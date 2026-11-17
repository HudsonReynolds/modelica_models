#!/usr/bin/env python3
"""Turn a picture of a gravitational lens into the matrix file that
GeneralRelativity.Examples.LensInversionFromFile reads.

The output is a MATLAB v4 file with one bordered matrix, lensImage[n + 1, n + 1],
the format that Modelica.Utilities.Streams.readRealMatrix understands without any
library (the v4 writer cannot hold several matrices):

    lensImage[2:, 1]   u, first image coordinate of the pixel centres, arcseconds
    lensImage[1, 2:]   v, second image coordinate, arcseconds (up)
    lensImage[2:, 2:]  surface brightness, background subtracted, peak scaled to 1
    lensImage[1, 1]    unused (0)

The brightness at (u[i], v[j]) is at row i + 1, column j + 1: the first index runs
along the picture's width, the second up the picture, as the Modelica functions expect.

Example:

    python tools/lens_image_to_mat.py arcs.png lens.mat --pixel-arcsec 0.05 --size 64 --crop 310 220 200 200

then in Modelica

    GeneralRelativity.Examples.LensInversionFromFile(imageFile = "lens.mat", zLens = 0.5, zSource = 2.0, psfArcsec = 0.08)

Only Pillow and NumPy are needed. The picture should show the arcs with the lens
galaxy's own light removed (or faint); a colour picture is converted to grey.
"""
import argparse
import struct
import sys

import numpy as np
from PIL import Image


def write_mat4(path, matrices):
    """Write {name: 2-D float array} as a little-endian MATLAB v4 file."""
    with open(path, "wb") as f:
        for name, m in matrices.items():
            m = np.ascontiguousarray(np.asarray(m, dtype="<f8"))
            if m.ndim != 2:
                raise ValueError(f"{name} must be a matrix")
            rows, cols = m.shape
            # type 0: little-endian IEEE double, full numeric matrix; imagf 0: real; name with its terminator.
            f.write(struct.pack("<iiiii", 0, rows, cols, 0, len(name) + 1))
            f.write(name.encode("ascii") + b"\0")
            f.write(m.tobytes(order="F"))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("picture", help="PNG, JPEG, TIFF or FITS-exported picture")
    p.add_argument("output", help="MATLAB v4 file to write (for example lens.mat)")
    p.add_argument("--pixel-arcsec", type=float, required=True, help="angular size of one picture pixel, arcseconds")
    p.add_argument("--size", type=int, default=64, help="resample to this many pixels per axis (default 64; the fit cost grows with the square)")
    p.add_argument("--crop", type=int, nargs=4, metavar=("LEFT", "TOP", "WIDTH", "HEIGHT"), help="crop box in picture pixels before resampling (default: the largest centred square)")
    p.add_argument("--background", choices=["median", "none"], default="median", help="subtract the median of the picture as the sky background (default) or nothing")
    p.add_argument("--invert", action="store_true", help="the picture is a negative (dark arcs on a light sky)")
    args = p.parse_args(argv)

    picture = Image.open(args.picture).convert("F")
    if args.crop:
        left, top, w, h = args.crop
        picture = picture.crop((left, top, left + w, top + h))
    else:
        w, h = picture.size
        s = min(w, h)
        picture = picture.crop(((w - s) // 2, (h - s) // 2, (w - s) // 2 + s, (h - s) // 2 + s))
    cropped_w, cropped_h = picture.size
    if cropped_w != cropped_h:
        sys.exit("the crop must be square so that the pixel scale is the same along both axes")
    # Box resampling averages the picture pixels inside each output pixel, which conserves surface brightness.
    picture = picture.resize((args.size, args.size), Image.BOX)
    arr = np.asarray(picture, dtype=float)
    if args.invert:
        arr = arr.max() - arr
    if args.background == "median":
        arr = arr - np.median(arr)
    arr = np.clip(arr, 0.0, None)
    peak = arr.max()
    if peak <= 0:
        sys.exit("the picture is empty after background subtraction")
    arr /= peak
    # Picture rows run from the top; image[i, j] wants i along the width and j upward.
    image = arr[::-1, :].T
    n = args.size
    scale = args.pixel_arcsec * cropped_w / n
    u = (np.arange(n) - (n - 1) / 2.0) * scale
    bordered = np.zeros((n + 1, n + 1))
    bordered[1:, 0] = u
    bordered[0, 1:] = u
    bordered[1:, 1:] = image
    write_mat4(args.output, {"lensImage": bordered})
    print(f"wrote {args.output}: {n} x {n} pixels of {scale:.4g} arcsec, field +/- {u[-1]:.3g} arcsec")
    print("Modelica: GeneralRelativity.Examples.LensInversionFromFile(imageFile = \"" + args.output.replace("\\", "/") + "\", zLens = ..., zSource = ..., psfArcsec = ...)")


if __name__ == "__main__":
    main()
