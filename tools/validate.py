#!/usr/bin/env python3
"""Comprueba el shader contra resultados exactos y contra el renderizador de referencia en CPU.

Uso (sin pantalla, desde la raíz del repo):
    cmake -S . -B build && cmake --build build --target reference_render
    xvfb-run -a python3 tools/validate.py --build build
"""
import argparse, math, os, subprocess, sys, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from render_shader import isco  # noqa: E402

failures = 0


def check(name, got, want, tol):
    global failures
    ok = abs(got - want) <= tol
    failures += not ok
    print(f"{'PASS' if ok else 'FAIL'} {name:<62} {got:12.6g}  esperado {want:10.6g}  (tol {tol:g})")


def shader_raw(tmp, cam, spin, size, target=(0, 0, 0), up=(0, 1, 0)):
    out = os.path.join(tmp, "shader.npy")
    subprocess.run([sys.executable, os.path.join(HERE, "render_shader.py"), "--cam", *map(str, cam),
                    "--target", *map(str, target), "--up", *map(str, up), "--size", *map(str, size),
                    "--set", f"spin={spin}", "--set", "viewMode=3", "--set", "turbulence=0",
                    "--raw", out, "--out", os.path.join(tmp, "shader.png")], check=True)
    return np.load(out)


def cpu_raw(build, tmp, cam, spin, size, target=(0, 0, 0), up=(0, 1, 0)):
    out = os.path.join(tmp, "cpu.f32")
    exe = os.path.join(build, "reference_render")
    subprocess.run([exe, str(size[0]), str(size[1]), *map(str, cam), str(spin), out, *map(str, target), *map(str, up)],
                   check=True)
    return np.fromfile(out, np.float32).reshape(size[1], size[0], 3)


def ut(r, a):
    r32 = r * math.sqrt(r)
    return (r32 + a) / (r ** 0.75 * math.sqrt(r32 - 3 * math.sqrt(r) + 2 * a))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default="build")
    args = ap.parse_args()
    tmp = tempfile.mkdtemp()

    print("== Vista desde el polo, cámara a r = 20M ==")
    for a in (0.0, 0.9):
        size = (801, 801)
        img = shader_raw(tmp, (0, 20, 0), a, size, up=(0, 0, -1))
        # Sombra: con a = 0 su radio angular es asin(3*sqrt(3)/r * sqrt(1 - 2/r))
        if a == 0.0:
            row = img[400, :, 2]
            dark = np.flatnonzero((row == 0) & (img[400, :, 3] == 0))
            r_px = (dark.max() - dark.min() + 1) / 2
            ang = math.degrees(math.atan(r_px / 400.5))
            want = math.degrees(math.asin(3 * math.sqrt(3) / 20 * math.sqrt(1 - 2 / 20)))
            check("K1 radio angular de la sombra, a=0 (grados)", ang, want, 0.15)
        # g de frente: 1 / (u^t E) con E = sqrt(1 - 2r/(r^2+a^2)) en el eje
        E = math.sqrt(1 - 2 * 20 / (400 + a * a))
        disk = img[..., 2] == 2
        g, r = img[..., 0][disk], img[..., 1][disk]
        want = np.array([1 / (ut(x, a) * E) for x in r])
        check(f"V2/K7 max error de g de frente, a={a}", float(np.max(np.abs(g / want - 1))), 0.0, 2e-3)
        for rr in (6.0, 10.0, 20.0):
            if rr > isco(a):
                print(f"     g esperado a {rr:4.1f}M: {1 / (ut(rr, a) * E):.4f}")
        # Temperatura: pico del flujo en 9.55M (a = 0)
        if a == 0.0:
            F = img[..., 3][disk]
            check("V4 radio del máximo de T, a=0", float(r[np.argmax(F)]), 9.548, 0.1)
            sel = np.abs(r - 15) < 0.05
            check("V4 T(15M)/Tmax, a=0", float(np.mean(F[sel]) ** 0.25), 0.884, 0.005)

    print("\n== K8: shader (float) frente a referencia en CPU (double), vista inclinada ==")
    for a, cam in ((0.0, (0, 2, 20)), (0.9, (0, 2, 20)), (0.9, (12, 6, -9)), (-0.7, (0, 8, 15))):
        size = (200, 160)
        s = shader_raw(tmp, cam, a, size)
        c = cpu_raw(args.build, tmp, cam, a, size)
        same = np.mean(s[..., 2] == c[..., 2])
        check(f"fracción de píxeles con el mismo resultado, a={a}, cam={cam}", float(same), 1.0, 0.01)
        both = (s[..., 2] == 2) & (c[..., 2] == 2)
        dg = np.abs(s[..., 0][both] / c[..., 0][both] - 1)
        check(f"percentil 99 del error relativo de g, a={a}, cam={cam}", float(np.percentile(dg, 99)), 0.0, 0.01)

    print(f"\n{'TODO OK' if failures == 0 else 'HAY FALLOS'} ({failures} fallos)")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
