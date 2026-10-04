// Física de Kerr compartida por las pruebas y el renderizador de referencia (doble precisión).
// simulator/blackhole.fs implementa exactamente lo mismo en GLSL; si cambias algo aquí, cámbialo allí.
//
// Unidades G = c = M = 1. Coordenadas Kerr-Schild cartesianas (X, Y, Z) con el spin a lo largo de +Z.
// En el simulador el eje "arriba" es y, y se usa (X, Y, Z) = (z, x, y), una permutación cíclica.
//
// Rayo de luz: hamiltoniano H = 1/2 g^{mu nu} p_mu p_nu = 0 con g^{mu nu} = eta^{mu nu} - f l^mu l^nu,
//   f = 2 r^3 / (r^4 + a^2 Z^2),  l_mu = (1, -(rX - aY)/(r^2+a^2), -(rY + aX)/(r^2+a^2), -Z/r).
// (Es la forma entrante con t -> -t y a -> -a.)
// p_t = -E es constante. Estado = posición X^i y momento covariante P_i.
#pragma once
#include <cmath>

namespace kerr {

struct Vec3 { double x, y, z; };
inline Vec3 operator+(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator-(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 operator*(Vec3 a, double s) { return {a.x * s, a.y * s, a.z * s}; }
inline double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(Vec3 a, Vec3 b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
inline double length(Vec3 a) { return std::sqrt(dot(a, a)); }
inline Vec3 normalize(Vec3 a) { return a * (1.0 / length(a)); }

inline double horizon(double a) { return 1.0 + std::sqrt(1.0 - a * a); }

// ISCO de Bardeen, Press y Teukolsky (1972). a > 0: disco co-rotante con el agujero.
inline double isco(double a) {
    double z1 = 1.0 + std::cbrt(1.0 - a * a) * (std::cbrt(1.0 + a) + std::cbrt(1.0 - a));
    double z2 = std::sqrt(3.0 * a * a + z1 * z1);
    double s = (a >= 0.0) ? 1.0 : -1.0;
    return 3.0 + z2 - s * std::sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2));
}

// Radio de Boyer-Lindquist (= radio Kerr-Schild) de un punto.
inline double radius(Vec3 X, double a) {
    double b = dot(X, X) - a * a;
    return std::sqrt(0.5 * (b + std::sqrt(b * b + 4.0 * a * a * X.z * X.z)));
}

struct Metric {   // f, l y sus gradientes en un punto
    double r, f;
    Vec3 l, gradR, gradF;
};

inline Metric metric(Vec3 X, double a) {
    Metric m;
    double r = radius(X, a), r2 = r * r, a2 = a * a;
    double q = r2 * r2 + a2 * X.z * X.z, s = r2 + a2;
    m.r = r;
    m.f = 2.0 * r2 * r / q;
    m.l = {-(r * X.x - a * X.y) / s, -(r * X.y + a * X.x) / s, -X.z / r};
    m.gradR = Vec3{X.x * r2 * r, X.y * r2 * r, X.z * r * s} * (1.0 / q);
    m.gradF = m.gradR * (2.0 * r2 * (3.0 * a2 * X.z * X.z - r2 * r2) / (q * q));
    m.gradF.z -= 4.0 * r2 * r * a2 * X.z / (q * q);
    return m;
}

// Gradiente de l·P respecto a la posición (P fijo).
inline Vec3 gradLdotP(Vec3 X, Vec3 P, double a, const Metric& m) {
    double r = m.r, s = r * r + a * a;
    double c = -(X.x * P.x + X.y * P.y) / s - 2.0 * r * (m.l.x * P.x + m.l.y * P.y) / s + X.z * P.z / (r * r);
    Vec3 g = m.gradR * c;
    g = g - Vec3{(r * P.x + a * P.y) / s, (r * P.y - a * P.x) / s, P.z / r};
    return g;
}

struct Ray { Vec3 X, P; double E; };   // E = -p_t

inline double hamiltonian(const Ray& y, double a) {
    Metric m = metric(y.X, a);
    double s = y.E + dot(m.l, y.P);
    return 0.5 * (-y.E * y.E + dot(y.P, y.P) - m.f * s * s);
}

// Derivadas respecto al parámetro afín: dX/dλ, dP/dλ y dt/dλ.
struct Deriv { Vec3 dX, dP; double dt; };
inline Deriv rhs(const Ray& y, double a) {
    Metric m = metric(y.X, a);
    double s = y.E + dot(m.l, y.P);
    Deriv d;
    d.dX = y.P - m.l * (m.f * s);
    d.dP = m.gradF * (0.5 * s * s) + gradLdotP(y.X, y.P, a, m) * (m.f * s);
    d.dt = y.E + m.f * s;
    return d;
}

// Un paso RK4 de tamaño h (negativo para trazar hacia el pasado). Devuelve dt acumulado.
inline double rk4(Ray& y, double a, double h) {
    auto add = [](const Ray& y, const Deriv& d, double k) { return Ray{y.X + d.dX * k, y.P + d.dP * k, y.E}; };
    Deriv k1 = rhs(y, a), k2 = rhs(add(y, k1, 0.5 * h), a), k3 = rhs(add(y, k2, 0.5 * h), a), k4 = rhs(add(y, k3, h), a);
    y.X = y.X + (k1.dX + k2.dX * 2.0 + k3.dX * 2.0 + k4.dX) * (h / 6.0);
    y.P = y.P + (k1.dP + k2.dP * 2.0 + k3.dP * 2.0 + k4.dP) * (h / 6.0);
    return (k1.dt + 2.0 * k2.dt + 2.0 * k3.dt + k4.dt) * (h / 6.0);
}

// Tamaño de paso del simulador: fino cerca del horizonte, proporcional a r lejos.
inline double stepSize(double r, double a, double scale = 1.0) {
    double rh = horizon(a);
    return scale * std::fmin(0.10 * (r - rh) + 0.02 * r, 0.25 * r);
}

// Producto escalar con la métrica g_{mu nu} = eta + f l l, vectores contravariantes (t, X).
struct V4 { double t; Vec3 s; };
inline double gdot(V4 u, V4 v, const Metric& m) {
    double lu = u.t + dot(m.l, u.s), lv = v.t + dot(m.l, v.s);
    return -u.t * v.t + dot(u.s, v.s) + m.f * lu * lv;
}

// Observador ZAMO (sin momento angular; estático si a = 0) y su base espacial ortonormal
// construida por Gram-Schmidt a partir de los ejes euclídeos de la cámara.
struct Tetrad { V4 u, fwd, up, right; };
inline Tetrad zamoTetrad(Vec3 X, double a, Vec3 fwd, Vec3 up, Vec3 right) {
    Metric m = metric(X, a);
    V4 xi{1.0, {0, 0, 0}}, psi{0.0, {-X.y, X.x, 0.0}};
    double pp = gdot(psi, psi, m);
    double w = pp > 1e-12 ? -gdot(xi, psi, m) / pp : 0.0;
    V4 u{1.0, psi.s * w};
    double nu = std::sqrt(-gdot(u, u, m));
    u = {u.t / nu, u.s * (1.0 / nu)};
    auto proj = [&](V4 v, V4 e, double sign) {   // quita la componente sobre e (sign = g(e,e))
        double c = gdot(v, e, m) / sign;
        return V4{v.t - c * e.t, v.s - e.s * c};
    };
    auto unit = [&](V4 v) { double n = std::sqrt(gdot(v, v, m)); return V4{v.t / n, v.s * (1.0 / n)}; };
    Tetrad T;
    T.u = u;
    T.fwd = unit(proj(V4{0, fwd}, u, -1));
    T.up = unit(proj(proj(V4{0, up}, u, -1), T.fwd, 1));
    T.right = unit(proj(proj(proj(V4{0, right}, u, -1), T.fwd, 1), T.up, 1));
    return T;
}

// Rayo que llega a la cámara desde la dirección local n = (adelante, arriba, derecha), energía local 1.
inline Ray cameraRay(Vec3 X, double a, const Tetrad& T, double nf, double nu, double nr) {
    Metric m = metric(X, a);
    V4 n{nf * T.fwd.t + nu * T.up.t + nr * T.right.t, T.fwd.s * nf + T.up.s * nu + T.right.s * nr};
    double norm = std::sqrt(nf * nf + nu * nu + nr * nr);
    V4 k{T.u.t - n.t / norm, T.u.s - n.s * (1.0 / norm)};    // el fotón viaja en dirección -n
    double lk = k.t + dot(m.l, k.s);
    // bajar índices: k_mu = eta_{mu nu} k^nu + f l_mu (l·k)
    double kt = -k.t + m.f * lk;
    Vec3 P = k.s + m.l * (m.f * lk);
    return Ray{X, P, -kt};
}

// Disco fino en órbitas circulares ecuatoriales (Boyer-Lindquist).
inline double diskOmega(double r, double a) { return 1.0 / (r * std::sqrt(r) + a); }
inline double diskUt(double r, double a) {
    double r32 = r * std::sqrt(r);
    return (r32 + a) / (std::pow(r, 0.75) * std::sqrt(r32 - 3.0 * std::sqrt(r) + 2.0 * a));
}
// Factor de corrimiento g = nu_obs / nu_emit para un punto del disco en (X, Y, 0).
inline double diskG(const Ray& y, double r, double a) {
    double pphi = y.X.x * y.P.y - y.X.y * y.P.x;
    return 1.0 / (diskUt(r, a) * (y.E - diskOmega(r, a) * pphi));
}

// Flujo de Page y Thorne (1974) con par nulo en la ISCO, normalización arbitraria.
inline double diskFlux(double r, double a) {
    double rin = isco(a);
    if (r <= rin) return 0.0;
    double x = std::sqrt(r), x0 = std::sqrt(rin), c = std::acos(a) / 3.0;
    double xs[3] = {2.0 * std::cos(c - M_PI / 3.0), 2.0 * std::cos(c + M_PI / 3.0), -2.0 * std::cos(c)};
    double s = x - x0 - 1.5 * a * std::log(x / x0);
    for (int i = 0; i < 3; i++) {
        double xi = xs[i], xj = xs[(i + 1) % 3], xk = xs[(i + 2) % 3];
        if (std::fabs(xi) < 1e-6) continue;   // a = 0: el término tiende a cero
        s -= 3.0 * (xi - a) * (xi - a) / (xi * (xi - xj) * (xi - xk)) * std::log((x - xi) / (x0 - xi));
    }
    return 1.5 * s / (x * x * x * x * (x * x * x - 3.0 * x + 2.0 * a));
}
inline double diskFluxMax(double a, double* rPeak = nullptr) {
    double rin = isco(a), best = 0.0, rb = rin;
    for (int i = 1; i <= 4000; i++) {
        double r = rin + i * 0.01 * (1.0 + 0.002 * i);
        double f = diskFlux(r, a);
        if (f > best) { best = f; rb = r; }
    }
    if (rPeak) *rPeak = rb;
    return best;
}

}  // namespace kerr
