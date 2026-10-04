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
    ap.add_argument("--size", type=int, nargs=2, default=[1000, 800])
    ap.add_argument("--time", type=float, default=0.0)
    ap.add_argument("--frames", type=int, default=1, help="accumulate N frames (iTime += 1/60 each)")
    ap.add_argument("--set", action="append", default=[], help="extra uniform, e.g. spin=0.9 or viewMode=1")
    ap.add_argument("--raw", help="also save float32 RGB as .npy")
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

    d, r, u = basis(a.cam, a.target)
    def setu(name, val):
        if name in prog:
            prog[name].value = val
    setu("camPos", tuple(a.cam)); setu("camDir", tuple(d)); setu("camRight", tuple(r)); setu("camUp", tuple(u))
    setu("resolution", (float(w), float(h)))
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
    acc = np.zeros((h, w, 3), np.float64)
    for i in range(a.frames):
        setu("iTime", a.time + i / 60.0)
        ctx.clear(0, 0, 0, 1)
        vao.render(moderngl.TRIANGLE_STRIP)
        img = np.frombuffer(fbo.read(components=4, dtype="f4"), np.float32).reshape(h, w, 4)
        acc += img[..., :3]
    acc /= a.frames
    acc = acc[::-1]  # GL origin is bottom-left
    if a.raw:
        np.save(a.raw, acc.astype(np.float32))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    Image.fromarray((np.clip(acc, 0, 1) * 255 + 0.5).astype(np.uint8)).save(a.out)

if __name__ == "__main__":
    main()
