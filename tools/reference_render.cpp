// Renderizador de referencia en CPU (doble precisión) con la misma física que simulator/blackhole.fs.
// Escribe, para cada píxel, (g, radio del impacto en el disco, tipo: 0 sombra / 1 cielo / 2 disco)
// en float32, igual que el modo de vista 3 del shader, para comparar los dos píxel a píxel.
//
// Uso: reference_render W H camX camY camZ spin salida.f32 [targetX targetY targetZ upX upY upZ] [escala_paso]
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "../common/KerrPhysics.h"
using namespace kerr;

static Vec3 toKS(Vec3 v) { return {v.z, v.x, v.y}; }

int main(int argc, char** argv) {
    if (argc < 8) {
        std::fprintf(stderr, "uso: %s W H camX camY camZ spin salida.f32 [target(3) up(3)] [escala]\n", argv[0]);
        return 1;
    }
    int W = std::atoi(argv[1]), H = std::atoi(argv[2]);
    Vec3 cam{std::atof(argv[3]), std::atof(argv[4]), std::atof(argv[5])};
    double a = std::atof(argv[6]);
    Vec3 target{0, 0, 0}, up{0, 1, 0};
    if (argc >= 14) {
        target = {std::atof(argv[8]), std::atof(argv[9]), std::atof(argv[10])};
        up = {std::atof(argv[11]), std::atof(argv[12]), std::atof(argv[13])};
    }
    double scale = argc >= 15 ? std::atof(argv[14]) : 1.0;
    const double fovY = M_PI / 2, rOut = 50.0;
    Vec3 d = normalize(target - cam), r = normalize(cross(d, up)), u = cross(r, d);
    Vec3 X0 = toKS(cam);
    Tetrad T = zamoTetrad(X0, a, toKS(d), toKS(u), toKS(r));
    double rh = horizon(a), rin = isco(a), rEsc = std::fmax(100.0, 2.0 * radius(X0, a));
    double tanHalf = std::tan(0.5 * fovY);
    int maxSteps = int(600 / scale) + 1;
    std::vector<float> out(size_t(W) * H * 3, 0.0f);
    for (int row = 0; row < H; row++) {
        for (int i = 0; i < W; i++) {
            double ux = (i + 0.5) / W * 2 - 1, uy = (H - 1 - row + 0.5) / H * 2 - 1;
            ux *= double(W) / H;
            Ray y = cameraRay(X0, a, T, 1.0, uy * tanHalf, ux * tanHalf);
            double pphi = y.X.x * y.P.y - y.X.y * y.P.x;
            float* px = &out[(size_t(row) * W + i) * 3];
            for (int s = 0; s < maxSteps; s++) {
                double rr = radius(y.X, a);
                if (rr < rh * 1.01) break;
                Deriv dd = rhs(y, a);
                if (rr > rEsc && dot(dd.dX, y.X) < 0) { px[2] = 1; break; }
                Vec3 Xold = y.X;
                rk4(y, a, -stepSize(rr, a, scale));
                if (Xold.z * y.X.z < 0) {
                    Vec3 Xh = Xold + (y.X - Xold) * (Xold.z / (Xold.z - y.X.z));
                    double rd = std::sqrt(std::fmax(Xh.x * Xh.x + Xh.y * Xh.y - a * a, 0.0));
                    if (rd > rin && rd < rOut) {
                        px[0] = float(1.0 / (diskUt(rd, a) * (y.E - diskOmega(rd, a) * pphi)));
                        px[1] = float(rd);
                        px[2] = 2;
                        break;
                    }
                }
            }
        }
    }
    FILE* f = std::fopen(argv[7], "wb");
    if (!f) { std::perror(argv[7]); return 1; }
    std::fwrite(out.data(), sizeof(float), out.size(), f);
    std::fclose(f);
    return 0;
}
