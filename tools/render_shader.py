#!/usr/bin/env python3
"""Render simulator/blackhole.fs offscreen to a PNG, with the same uniforms main.cpp sends.

Needs: pip install moderngl numpy pillow. On a machine without a display, run it under xvfb-run.
Example: python3 tools/render_shader.py --cam 0 2 20 --out out/default.png
"""
import argparse, math, os, sys
import numpy as np
import moderngl
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
VERT = """#version 330
in vec2 in_pos;
void main() { gl_Position = vec4(in_pos, 0.0, 1.0); }
"""

def isco(a):
    z1 = 1 + (1 - a * a) ** (1 / 3) * ((1 + a) ** (1 / 3) + (1 - a) ** (1 / 3))
    z2 = math.sqrt(3 * a * a + z1 * z1)
    return 3 + z2 - (1 if a >= 0 else -1) * math.sqrt(max((3 - z1) * (3 + z1 + 2 * z2), 0))

def disk_flux(r, a):
    """Page-Thorne, same as diskFlux in common/KerrPhysics.h."""
    rin = isco(a)
    if r <= rin:
        return 0.0
    x, x0, c = math.sqrt(r), math.sqrt(rin), math.acos(a) / 3
    xs = [2 * math.cos(c - math.pi / 3), 2 * math.cos(c + math.pi / 3), -2 * math.cos(c)]
    s = x - x0 - 1.5 * a * math.log(x / x0)
    for i in range(3):
        xi, xj, xk = xs[i], xs[(i + 1) % 3], xs[(i + 2) % 3]
        if abs(xi) < 1e-6:
            continue
        s -= 3 * (xi - a) ** 2 / (xi * (xi - xj) * (xi - xk)) * math.log((x - xi) / (x0 - xi))
    return 1.5 * s / (x ** 4 * (x ** 3 - 3 * x + 2 * a))

def disk_flux_max(a):
    rin = isco(a)
    return max(disk_flux(rin + i * 0.01 * (1 + 0.002 * i), a) for i in range(1, 4001))

DEFAULTS = {"spin": 0.0, "fovY": math.pi / 2, "Tmax": 6500.0, "exposure": 0.6, "viewMode": 0,
            "turbulence": 1, "timeScale": 10.0}

def basis(cam, target, up=(0.0, 1.0, 0.0)):
    cam, target, up = map(np.asarray, (cam, target, up))
    d = target - cam; d = d / np.linalg.norm(d)
    r = np.cross(d, up); r = r / np.linalg.norm(r)
    u = np.cross(r, d)
    return d, r, u

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shader", default=os.path.join(HERE, "..", "simulator", "blackhole.fs"))
    ap.add_argument("--cam", type=float, nargs=3, default=[0.0, 2.0, 20.0])
    ap.add_argument("--target", type=float, nargs=3, default=[0.0, 0.0, 0.0])
    ap.add_argument("--up", type=float, nargs=3, default=[0.0, 1.0, 0.0])
    ap.add_argument("--size", type=int, nargs=2, default=[1000, 800])
    ap.add_argument("--time", type=float, default=0.0)
    ap.add_argument("--frames", type=int, default=1, help="accumulate N frames (iTime += 1/60 each)")
    ap.add_argument("--set", action="append", default=[], help="extra uniform, e.g. spin=0.9 or viewMode=1")
    ap.add_argument("--raw", help="also save float32 RGBA as .npy (viewMode=3 gives g, r, type, F)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    ctx = moderngl.create_standalone_context(require=330)
    src = open(a.shader).read()
    prog = ctx.program(vertex_shader=VERT, fragment_shader=src)
    quad = ctx.buffer(np.array([-1, -1, 1, -1, -1, 1, 1, 1], dtype="f4").tobytes())
    vao = ctx.vertex_array(prog, [(quad, "2f", "in_pos")])
    w, h = a.size
    fbo = ctx.framebuffer(color_attachments=[ctx.texture((w, h), 4, dtype="f4")])
    fbo.use()

    d, r, u = basis(a.cam, a.target, a.up)
    def setu(name, val):
        if name in prog:
            prog[name].value = val
    setu("camPos", tuple(a.cam)); setu("camDir", tuple(d)); setu("camRight", tuple(r)); setu("camUp", tuple(u))
    setu("resolution", (float(w), float(h)))
    for k, v in DEFAULTS.items():
        if k in prog:
            prog[k].value = v
    for kv in a.set:
        k, v = kv.split("=")
        vals = [float(x) for x in v.split(",")]
        if k not in prog:
            print(f"warning: uniform {k} not in shader", file=sys.stderr); continue
        cur = prog[k]
        if cur.dimension == 1:
            cur.value = int(vals[0]) if cur.fmt.endswith("i") else vals[0]
        else:
            cur.value = tuple(vals)
    if "fluxMax" in prog:
        prog["fluxMax"].value = disk_flux_max(prog["spin"].value if "spin" in prog else 0.0)
    acc = np.zeros((h, w, 4), np.float64)
    for i in range(a.frames):
        setu("iTime", a.time + i / 60.0)
        ctx.clear(0, 0, 0, 1)
        vao.render(moderngl.TRIANGLE_STRIP)
        img = np.frombuffer(fbo.read(components=4, dtype="f4"), np.float32).reshape(h, w, 4)
        acc += img
    acc /= a.frames
    acc = acc[::-1]  # GL origin is bottom-left
    if a.raw:
        np.save(a.raw, acc.astype(np.float32))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    Image.fromarray((np.clip(acc[..., :3], 0, 1) * 255 + 0.5).astype(np.uint8)).save(a.out)

if __name__ == "__main__":
    main()
