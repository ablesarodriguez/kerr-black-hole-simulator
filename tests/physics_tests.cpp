// Pruebas numéricas de common/KerrPhysics.h contra resultados analíticos exactos.
// Compilar: cmake --build build --target physics_tests && ./build/physics_tests
#include <cstdio>
#include <cmath>
#include <initializer_list>
#include "../common/KerrPhysics.h"
using namespace kerr;

static int failures = 0;
static void check(const char* name, double got, double want, double tolRel) {
    double err = std::fabs(got - want) / std::fmax(std::fabs(want), 1e-12);
    bool ok = err <= tolRel;
    if (!ok) failures++;
    std::printf("%s %-58s got %.6f  want %.6f  (err %.2e)\n", ok ? "PASS" : "FAIL", name, got, want, err);
}

enum Outcome { CAPTURED, ESCAPED, DISK, MAXED };
struct Trace { Outcome out; int steps; double diskR; double maxDriftH, maxDriftL; };

// Traza hacia el pasado como el shader. stopAtDisk: parar en el primer cruce de Z = 0 fuera de la ISCO.
static Trace trace(Ray y, double a, double scale, bool stopAtDisk, double rEsc = 1100.0, int maxSteps = 200000) {
    Trace t{MAXED, 0, 0, 0, 0};
    double L0 = y.X.x * y.P.y - y.X.y * y.P.x, rh = horizon(a), rin = isco(a);
    for (; t.steps < maxSteps; t.steps++) {
        double r = radius(y.X, a);
        if (r < rh * 1.01) { t.out = CAPTURED; return t; }
        Deriv d = rhs(y, a);
        if (r > rEsc && dot(d.dX, y.X) < 0) { t.out = ESCAPED; return t; }   // hacia atrás se aleja
        Vec3 Xold = y.X;
        rk4(y, a, -stepSize(r, a, scale));
        t.maxDriftH = std::fmax(t.maxDriftH, std::fabs(hamiltonian(y, a)));
        t.maxDriftL = std::fmax(t.maxDriftL, std::fabs((y.X.x * y.P.y - y.X.y * y.P.x) - L0));
        if (stopAtDisk && Xold.z * y.X.z <= 0.0) {
            double w = Xold.z / (Xold.z - y.X.z);
            Vec3 Xh = Xold + (y.X - Xold) * w;
            double rr = std::sqrt(std::fmax(Xh.x * Xh.x + Xh.y * Xh.y - a * a, 0.0));
            if (rr > rin) { t.out = DISK; t.diskR = rr; return t; }
        }
    }
    return t;
}

// Cámara ZAMO en el ecuador mirando al agujero; rayo con ángulo alpha (rad) hacia "derecha".
static Ray equatorialRay(double rc, double a, double alpha) {
    Vec3 X{rc, 0, 0}, fwd{-1, 0, 0}, up{0, 0, 1}, right{0, 1, 0};
    Tetrad T = zamoTetrad(X, a, fwd, up, right);
    return cameraRay(X, a, T, std::cos(alpha), 0.0, std::sin(alpha));
}
static double criticalAngle(double rc, double a, double lo, double hi, double scale) {   // lo capturado
    for (int i = 0; i < 60; i++) {
        double mid = 0.5 * (lo + hi);
        (trace(equatorialRay(rc, a, mid), a, scale, false).out == CAPTURED ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}
// Bardeen (1973): órbitas circulares ecuatoriales de fotones r(r-3)^2 = 4a^2, xi = L/E.
static double bardeenEdge(double a, bool prograde) {
    if (a == 0.0) return prograde ? 3.0 * std::sqrt(3.0) : -3.0 * std::sqrt(3.0);
    double lo = prograde ? 1.0 : 3.0, hi = prograde ? 3.0 : 4.0;
    auto g = [&](double r) { return r * (r - 3) * (r - 3) - 4 * a * a; };
    for (int i = 0; i < 200; i++) { double m = 0.5 * (lo + hi); (g(lo) * g(m) <= 0 ? hi : lo) = m; }
    double r = 0.5 * (lo + hi);
    return -(r * r * r - 3 * r * r + a * a * r + a * a) / (a * (r - 1));
}

int main() {
    std::printf("== Gradiente analítico frente a diferencias finitas ==\n");
    {
        double worst = 0;
        for (double a : {0.0, 0.6, 0.95, -0.8})
            for (int k = 0; k < 50; k++) {
                Vec3 X{std::sin(k * 1.7) * (3 + k % 7), std::cos(k * 2.3) * (2 + k % 5), std::sin(k * 0.9) * (1 + k % 4)};
                Ray y{X, {std::sin(k * 1.1), std::cos(k * 0.7), std::sin(k * 0.3 + 1)}, 1.0};
                Deriv d = rhs(y, a);
                double e = 1e-6;
                for (int i = 0; i < 3; i++) {
                    Ray yp = y, ym = y;
                    (&yp.X.x)[i] += e; (&ym.X.x)[i] -= e;
                    double num = -(hamiltonian(yp, a) - hamiltonian(ym, a)) / (2 * e);
                    worst = std::fmax(worst, std::fabs(num - (&d.dP.x)[i]) / (1 + std::fabs(num)));
                    Ray qp = y, qm = y;
                    (&qp.P.x)[i] += e; (&qm.P.x)[i] -= e;
                    num = (hamiltonian(qp, a) - hamiltonian(qm, a)) / (2 * e);
                    worst = std::fmax(worst, std::fabs(num - (&d.dX.x)[i]) / (1 + std::fabs(num)));
                }
            }
        check("1 + max error relativo de dX/dlambda y dP/dlambda", 1.0 + worst, 1.0, 1e-7);
    }

    std::printf("\n== Radios característicos ==\n");
    check("ISCO a=0", isco(0), 6.0, 1e-12);
    check("ISCO a=0.9 progrado", isco(0.9), 2.320883, 1e-6);
    check("ISCO a=-0.9 (retrógrado)", isco(-0.9), 8.717352, 1e-6);
    check("horizonte a=0.9", horizon(0.9), 1.435890, 1e-6);
    double rp;
    diskFluxMax(0.0, &rp);
    check("pico del flujo de Page-Thorne a=0", rp, 9.548, 2e-3);
    check("T(15M)/Tmax a=0", std::pow(diskFlux(15, 0) / diskFluxMax(0), 0.25), 0.884, 2e-3);

    std::printf("\n== Observador: el rayo sale nulo y con energía local 1 ==\n");
    for (double a : {0.0, 0.9}) {
        Vec3 X{3.0, 4.0, 5.0};
        Tetrad T = zamoTetrad(X, a, normalize(Vec3{-3, -4, -5}), {0, 0, 1}, normalize(cross(Vec3{-3, -4, -5}, {0, 0, 1})));
        Ray y = cameraRay(X, a, T, 0.3, -0.5, 0.8);
        char buf[80]; std::snprintf(buf, sizeof buf, "|H| del rayo inicial, a=%.1f", a);
        check(buf, std::fabs(hamiltonian(y, a)) + 1.0, 1.0, 1e-12);
    }

    std::printf("\n== K1: sombra vista por un observador estático en r = 20M, a = 0 ==\n");
    {
        double alpha = criticalAngle(20.0, 0.0, 0.0, 0.6, 1.0);
        double want = std::asin(3.0 * std::sqrt(3.0) / 20.0 * std::sqrt(1.0 - 2.0 / 20.0));
        check("radio angular de la sombra (grados)", alpha * 180 / M_PI, want * 180 / M_PI, 1e-3);
    }

    std::printf("\n== K2: bordes de la sombra de Kerr (Bardeen 1973), cámara en r = 1000M ==\n");
    for (double scale : {0.1, 1.0})
        for (double a : {0.0, 0.5, 0.9, 0.95}) {
            for (int side = 0; side < 2; side++) {
                double s = side == 0 ? 1.0 : -1.0;
                double al = criticalAngle(1000.0, a, 0.0, s * 0.02, scale);
                Ray y = equatorialRay(1000.0, a, al);
                double b = (y.X.x * y.P.y - y.X.y * y.P.x) / y.E;
                // el rayo trazado hacia atrás con alpha > 0 es un fotón retrógrado (L < 0) visto desde +X
                double want = bardeenEdge(a, b > 0);
                char buf[96];
                std::snprintf(buf, sizeof buf, "a=%.2f borde %s, paso x%.1f", a, b > 0 ? "progrado " : "retrógrado", scale);
                check(buf, b, want, scale < 1 ? 1e-4 : 2e-3);
            }
        }

    std::printf("\n== K5: conservación a lo largo de un rayo rasante (a=0.9, paso del simulador) ==\n");
    {
        double al = criticalAngle(20.0, 0.9, 0.0, 0.6, 1.0);
        Trace t = trace(equatorialRay(20.0, 0.9, al * 1.01), 0.9, 1.0, false, 100.0);
        std::printf("     pasos: %d\n", t.steps);
        check("deriva max |H|", t.maxDriftH + 1.0, 1.0, 1e-4);
        check("deriva max de L_z", t.maxDriftL + 1.0, 1.0, 1e-4);
    }

    std::printf("\n== K4: a = 0 en Kerr-Schild coincide con la fuerza efectiva (radio del disco) ==\n");
    {
        double worst = 0;
        int n = 0, steps = 0;
        for (int i = 0; i < 40; i++) {
            Vec3 X{0, 20.0 * std::sin(0.1), 20.0 * std::cos(0.1)};   // cámara del simulador (0,2,20) aprox., ejes KS
            Vec3 fwd = normalize(X * -1.0), up{0, 0, 1};
            Vec3 right = normalize(cross(fwd, up)); up = cross(right, fwd);
            Tetrad T = zamoTetrad(X, 0.0, fwd, up, right);
            double nu = -0.35 + 0.7 * (i % 8) / 7.0, nr = -0.35 + 0.7 * (i / 8) / 4.0;
            Ray y = cameraRay(X, 0.0, T, 1.0, nu, nr);
            Trace t = trace(y, 0.0, 1.0, true, 100.0);
            if (t.out != DISK) continue;
            steps += t.steps;
            // referencia: fuerza efectiva x'' = -3 h^2 x / r^5 con dirección euclídea equivalente
            Deriv d = rhs(y, 0.0);
            Vec3 v = normalize(d.dX * -1.0), x = X;
            double h = 1e-3;
            for (int k = 0; k < 2000000; k++) {
                auto acc = [](Vec3 p, Vec3 v) { Vec3 hh = cross(p, v); double r = length(p); return p * (-3.0 * dot(hh, hh) / std::pow(r, 5)); };
                Vec3 xo = x;
                Vec3 k1v = acc(x, v), k1x = v, k2v = acc(x + k1x * (h / 2), v + k1v * (h / 2)), k2x = v + k1v * (h / 2);
                Vec3 k3v = acc(x + k2x * (h / 2), v + k2v * (h / 2)), k3x = v + k2v * (h / 2);
                Vec3 k4v = acc(x + k3x * h, v + k3v * h), k4x = v + k3v * h;
                v = normalize(v + (k1v + k2v * 2 + k3v * 2 + k4v) * (h / 6));
                x = x + (k1x + k2x * 2 + k3x * 2 + k4x) * (h / 6);
                if (xo.z * x.z <= 0) {
                    Vec3 xh = xo + (x - xo) * (xo.z / (xo.z - x.z));
                    double rr = std::sqrt(xh.x * xh.x + xh.y * xh.y);
                    if (rr > 6.0) { worst = std::fmax(worst, std::fabs(rr - t.diskR) / rr); n++; break; }
                }
            }
        }
        std::printf("     %d rayos al disco, %.0f pasos de media\n", n, n ? double(steps) / n : 0.0);
        check("max error relativo del radio de impacto en el disco", worst + 1.0, 1.0, 2e-3);
    }

    std::printf("\n%s (%d fallos)\n", failures ? "HAY FALLOS" : "TODO OK", failures);
    return failures ? 1 : 0;
}
