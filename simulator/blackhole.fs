#version 330
// Agujero negro de Kerr con disco de acreción fino, trazado por geodésicas nulas.
//
// Unidades G = c = M = 1 (la cámara se mide en M). La física es la de common/KerrPhysics.h,
// que tiene pruebas numéricas en tests/physics_tests.cpp; si cambias algo aquí, cámbialo allí.
//
//  - Rayos: hamiltoniano en coordenadas Kerr-Schild salientes cartesianas, RK4, trazado hacia el pasado.
//  - Cámara: observador ZAMO (estático si spin = 0) con base ortonormal local.
//  - Disco: fino en el plano y = 0, desde la ISCO, flujo de Page-Thorne, cuerpo negro a g·T.
//  - Cielo: estrellas de cuerpo negro desplazadas al azul por la gravedad en la cámara.

uniform vec3 camPos;
uniform vec2 resolution;
uniform vec3 camDir;
uniform vec3 camRight;
uniform vec3 camUp;
uniform float iTime;

uniform float spin;        // a/M, entre -0.998 y 0.998 (negativo: disco retrógrado)
uniform float fovY;        // campo de visión vertical en radianes
uniform float Tmax;        // temperatura máxima del disco en K (en el pico del flujo)
uniform float fluxMax;     // máximo del flujo de Page-Thorne para este spin (lo calcula main.cpp)
uniform float exposure;
uniform int viewMode;      // 0 imagen, 1 mapa de g, 2 mapa de T emitida, 3 datos crudos (g, r, tipo) para validar
uniform bool turbulence;
uniform float timeScale;   // M de tiempo coordenado por segundo de reloj

out vec4 finalColor;

const float PI = 3.14159265;
const float R_DISK_OUT = 25.0;   // el disco se vuelve transparente entre 25M y 50M
const float R_DISK_FADE = 50.0;
const int MAX_STEPS = 600;

// ---------- Ruido ----------
float hash(vec2 p) {
    return fract(sin(dot(p, vec2(127.1, 311.7))) * 43758.5453123);
}

float noise(vec2 p) {
    vec2 i = floor(p);
    vec2 f = fract(p);
    float a = hash(i);
    float b = hash(i + vec2(1.0, 0.0));
    float c = hash(i + vec2(0.0, 1.0));
    float d = hash(i + vec2(1.0, 1.0));
    vec2 u = f*f*(3.0-2.0*f);
    return mix(a, b, u.x) + (c-a)*u.y*(1.0-u.x) + (d-b)*u.x*u.y;
}

float fbm(vec2 p) {
    float value = 0.0;
    float amp = 0.5;
    for (int i = 0; i < 3; i++) {
        value += amp * noise(p);
        p *= 2.02;
        amp *= 0.5;
    }
    return value;
}

// ---------- Cuerpo negro ----------
// Funciones de color CIE 1931 (ajuste de Wyman, Sloan y Shirley 2013).
float lobe(float x, float mu, float s1, float s2) {
    float t = (x - mu) / (x < mu ? s1 : s2);
    return exp(-0.5 * t * t);
}
vec3 cieXYZ(float nm) {
    return vec3(1.056 * lobe(nm, 599.8, 37.9, 31.0) + 0.362 * lobe(nm, 442.0, 16.0, 26.7) - 0.065 * lobe(nm, 501.1, 20.4, 26.2),
                0.821 * lobe(nm, 568.8, 46.9, 40.5) + 0.286 * lobe(nm, 530.9, 16.3, 31.1),
                1.217 * lobe(nm, 437.0, 11.8, 36.0) + 0.681 * lobe(nm, 459.0, 26.0, 13.8));
}
// Radiancia de cuerpo negro a temperatura T integrada con las funciones CIE (unidades relativas).
vec3 blackbodyXYZ(float T) {
    vec3 xyz = vec3(0.0);
    for (int i = 0; i < 41; i++) {
        float nm = 380.0 + 10.0 * float(i);
        float um = nm * 1.0e-3;
        float planck = 1.0 / (um * um * um * um * um * (exp(min(14387.77 / (um * T), 80.0)) - 1.0));
        xyz += planck * cieXYZ(nm);
    }
    return xyz;
}
vec3 xyzToLinearSRGB(vec3 c) {
    return vec3( 3.2406 * c.x - 1.5372 * c.y - 0.4986 * c.z,
                -0.9689 * c.x + 1.8758 * c.y + 0.0415 * c.z,
                 0.0557 * c.x - 0.2040 * c.y + 1.0570 * c.z);
}
// Color observado de un cuerpo negro emitido a T y visto con corrimiento g, relativo al brillo
// de un cuerpo negro a Tref visto sin corrimiento. I_nu / nu^3 se conserva, así que el espectro
// observado es exactamente un cuerpo negro a g·T: de ahí salen a la vez el color y el brillo.
vec3 observedBlackbody(float T, float g, float Tref) {
    return max(xyzToLinearSRGB(blackbodyXYZ(g * T) / blackbodyXYZ(Tref).y), vec3(0.0));
}

// ---------- Métrica de Kerr en Kerr-Schild salientes (spin a lo largo de +Z) ----------
// Coordenadas del shader (x, y, z) con y hacia arriba -> Kerr-Schild (X, Y, Z) = (z, x, y).
vec3 toKS(vec3 v)   { return vec3(v.z, v.x, v.y); }
vec3 fromKS(vec3 v) { return vec3(v.y, v.z, v.x); }

float horizonRadius(float a) { return 1.0 + sqrt(1.0 - a * a); }

float iscoRadius(float a) {
    float z1 = 1.0 + pow(1.0 - a * a, 1.0 / 3.0) * (pow(1.0 + a, 1.0 / 3.0) + pow(1.0 - a, 1.0 / 3.0));
    float z2 = sqrt(3.0 * a * a + z1 * z1);
    float s = a >= 0.0 ? 1.0 : -1.0;
    return 3.0 + z2 - s * sqrt(max((3.0 - z1) * (3.0 + z1 + 2.0 * z2), 0.0));
}

float radiusKS(vec3 X, float a) {
    float b = dot(X, X) - a * a;
    return sqrt(0.5 * (b + sqrt(b * b + 4.0 * a * a * X.z * X.z)));
}

struct Metric { float r; float f; vec3 l; vec3 gradR; vec3 gradF; };

Metric metricAt(vec3 X, float a) {
    Metric m;
    float r = radiusKS(X, a);
    float r2 = r * r, a2 = a * a;
    float q = r2 * r2 + a2 * X.z * X.z, s = r2 + a2;
    m.r = r;
    m.f = 2.0 * r2 * r / q;
    m.l = vec3(-(r * X.x - a * X.y) / s, -(r * X.y + a * X.x) / s, -X.z / r);
    m.gradR = vec3(X.x * r2 * r, X.y * r2 * r, X.z * r * s) / q;
    m.gradF = m.gradR * (2.0 * r2 * (3.0 * a2 * X.z * X.z - r2 * r2) / (q * q));
    m.gradF.z -= 4.0 * r2 * r * a2 * X.z / (q * q);
    return m;
}

vec3 gradLdotP(vec3 X, vec3 P, float a, Metric m) {
    float r = m.r, s = r * r + a * a;
    float c = -(X.x * P.x + X.y * P.y) / s - 2.0 * r * (m.l.x * P.x + m.l.y * P.y) / s + X.z * P.z / (r * r);
    return m.gradR * c - vec3((r * P.x + a * P.y) / s, (r * P.y - a * P.x) / s, P.z / r);
}

// dX/dλ, dP/dλ y dt/dλ para el hamiltoniano H = 1/2 [-E^2 + |P|^2 - f (E + l·P)^2].
void geodesic(vec3 X, vec3 P, float E, float a, out vec3 dX, out vec3 dP, out float dt) {
    Metric m = metricAt(X, a);
    float s = E + dot(m.l, P);
    dX = P - m.l * (m.f * s);
    dP = m.gradF * (0.5 * s * s) + gradLdotP(X, P, a, m) * (m.f * s);
    dt = E + m.f * s;
}

float stepSize(float r, float rh) {
    return min(0.10 * (r - rh) + 0.02 * r, 0.25 * r);
}

// Producto escalar con g_{mu nu} = eta + f l l para vectores contravariantes (t, X).
float gdot(vec4 u, vec4 v, Metric m) {
    float lu = u.x + dot(m.l, u.yzw), lv = v.x + dot(m.l, v.yzw);
    return -u.x * v.x + dot(u.yzw, v.yzw) + m.f * lu * lv;
}
// Gram-Schmidt con la métrica: quita a v sus componentes sobre u (temporal) y sobre e1, e2.
vec4 orthonormalize(vec4 v, vec4 u, vec4 e1, vec4 e2, int n, Metric m) {
    v += gdot(v, u, m) * u;              // g(u,u) = -1
    if (n > 0) v -= gdot(v, e1, m) * e1;
    if (n > 1) v -= gdot(v, e2, m) * e2;
    return v / sqrt(gdot(v, v, m));
}

// ---------- Disco ----------
float diskOmega(float r, float a) { return 1.0 / (r * sqrt(r) + a); }
float diskUt(float r, float a) {
    float r32 = r * sqrt(r);
    return (r32 + a) / (pow(r, 0.75) * sqrt(r32 - 3.0 * sqrt(r) + 2.0 * a));
}
// Flujo de Page y Thorne (1974), par nulo en la ISCO (normalización arbitraria).
float diskFlux(float r, float a, float rin) {
    if (r <= rin) return 0.0;
    float x = sqrt(r), x0 = sqrt(rin), c = acos(clamp(a, -1.0, 1.0)) / 3.0;
    float xs[3] = float[3](2.0 * cos(c - PI / 3.0), 2.0 * cos(c + PI / 3.0), -2.0 * cos(c));
    float s = x - x0 - 1.5 * a * log(x / x0);
    for (int i = 0; i < 3; i++) {
        float xi = xs[i], xj = xs[(i + 1) % 3], xk = xs[(i + 2) % 3];
        if (abs(xi) < 1.0e-4) continue;   // a = 0: el término tiende a cero
        s -= 3.0 * (xi - a) * (xi - a) / (xi * (xi - xj) * (xi - xk)) * log((x - xi) / (x0 - xi));
    }
    return max(1.5 * s / (x * x * x * x * (x * x * x - 3.0 * x + 2.0 * a)), 0.0);
}

// Turbulencia arrastrada por la rotación diferencial. Dos copias desfasadas medio periodo se
// mezclan para que el patrón no se enrolle indefinidamente.
float diskTurbulence(float r, float phi, float omega, float t) {
    const float PERIOD = 300.0;
    float tp = t / PERIOD, n = 0.0;
    for (int k = 0; k < 2; k++) {
        float ph = fract(tp + 0.5 * float(k));
        float w = 1.0 - abs(2.0 * ph - 1.0);
        float seed = floor(tp + 0.5 * float(k)) + float(k) * 0.37;
        float ang = phi - omega * PERIOD * ph;
        vec2 p = r * vec2(cos(ang), sin(ang)) * 0.18 + vec2(seed * 17.3, seed * 5.1);
        n += w * fbm(p);
    }
    return 0.4 + 1.2 * n;
}

// ---------- Cielo ----------
vec3 renderSky(vec3 dir, float gSky) {
    vec2 sph = vec2(atan(dir.z, dir.x), acos(clamp(dir.y, -1.0, 1.0)));
    vec2 grid = sph * vec2(40.0, 40.0);
    vec2 id = floor(grid);
    vec2 f = fract(grid) - 0.5;

    float h = hash(id);
    vec3 star = vec3(0.0);
    if (h > 0.985) {
        vec2 jitter = vec2(hash(id + 3.7), hash(id + 91.3)) - 0.5;
        float d = length(f - jitter * 0.6);
        float brightness = (h - 0.985) / 0.015;
        float shape = smoothstep(0.12, 0.0, d) * brightness;
        if (shape > 0.0) {
            float Tstar = 3000.0 + 22000.0 * pow(hash(id + 7.1), 2.0);
            star = 0.6 * shape * observedBlackbody(Tstar, gSky, Tstar);
        }
    }
    float band = pow(max(0.0, 1.0 - abs(dir.y)), 6.0);
    vec3 bg = vec3(0.0002, 0.00025, 0.0005) + band * vec3(0.004, 0.0037, 0.0047);
    // El fondo difuso se trata como gris: su brillo escala como g^4.
    return bg * pow(gSky, 4.0) + star;
}

// ---------- Mapas de depuración ----------
vec3 divergingMap(float g) {   // g < 1 rojo, g = 1 blanco, g > 1 azul
    float t = clamp((g - 1.0) / 0.7, -1.0, 1.0);
    return t < 0.0 ? mix(vec3(1.0), vec3(0.8, 0.1, 0.05), -t) : mix(vec3(1.0), vec3(0.1, 0.25, 0.9), t);
}
vec3 heatMap(float x) {
    x = clamp(x, 0.0, 1.0);
    return clamp(vec3(3.0 * x, 3.0 * x - 1.0, 3.0 * x - 2.0), 0.0, 1.0);
}

vec3 tonemap(vec3 x) {   // ACES aproximado (Narkowicz) + codificación sRGB
    x = clamp((x * (2.51 * x + 0.03)) / (x * (2.43 * x + 0.59) + 0.14), 0.0, 1.0);
    return pow(x, vec3(1.0 / 2.2));
}

void main() {
    float a = clamp(spin, -0.998, 0.998);
    float rh = horizonRadius(a);
    float rin = iscoRadius(a);

    // Dirección del píxel en el marco local de la cámara
    vec2 uv = gl_FragCoord.xy / resolution.xy;
    vec2 jitter = vec2(hash(uv + iTime), hash(uv + iTime + 1.0)) - 0.5;
    uv = uv * 2.0 - 1.0 + (viewMode == 3 ? vec2(0.0) : jitter * 0.0005);
    uv.x *= resolution.x / resolution.y;
    float tanHalf = tan(0.5 * fovY);
    vec3 n = normalize(vec3(1.0, uv.y * tanHalf, uv.x * tanHalf));   // (adelante, arriba, derecha)

    // Base ortonormal del observador ZAMO en la posición de la cámara
    vec3 X = toKS(camPos);
    Metric m0 = metricAt(X, a);
    vec4 xi = vec4(1.0, 0.0, 0.0, 0.0), psi = vec4(0.0, -X.y, X.x, 0.0);
    float pp = gdot(psi, psi, m0);
    float w = pp > 1.0e-8 ? -gdot(xi, psi, m0) / pp : 0.0;
    vec4 u = xi + w * psi;
    float uu = -gdot(u, u, m0);
    if (uu <= 0.0 || m0.r < rh) { finalColor = vec4(0.0, 0.0, 0.0, 1.0); return; }   // cámara dentro del horizonte
    u /= sqrt(uu);
    vec4 eF = orthonormalize(vec4(0.0, toKS(camDir)), u, u, u, 0, m0);
    vec4 eU = orthonormalize(vec4(0.0, toKS(camUp)), u, eF, eF, 1, m0);
    vec4 eR = orthonormalize(vec4(0.0, toKS(camRight)), u, eF, eU, 2, m0);

    // Fotón que llega a la cámara viajando en la dirección -n, con energía local 1
    vec4 k = u - (n.x * eF + n.y * eU + n.z * eR);
    float lk = k.x + dot(m0.l, k.yzw);
    float E = k.x - m0.f * lk;                // E = -p_t
    vec3 P = k.yzw + m0.l * (m0.f * lk);      // momento covariante espacial
    float pphi = X.x * P.y - X.y * P.x;       // L_z, constante del movimiento

    float rEsc = max(100.0, 2.0 * m0.r);
    float tObs = iTime * timeScale;
    float tAcc = 0.0;                         // tiempo coordenado recorrido hacia el pasado (<= 0)

    vec3 color = vec3(0.0);
    float transmit = 1.0;
    bool escaped = false;
    bool haveDebug = false;
    vec3 debug = vec3(0.0);
    vec4 raw = vec4(0.0);                     // modo 3: (g, radio del impacto, 0 sombra / 1 cielo / 2 disco, F/Fmax)
    vec3 skyDir = vec3(0.0);

    for (int i = 0; i < MAX_STEPS; i++) {
        vec3 k1x, k1p; float k1t;
        geodesic(X, P, E, a, k1x, k1p, k1t);
        float r = radiusKS(X, a);
        if (r < rh * 1.01) break;                                  // viene del horizonte: sombra
        if (r > rEsc && dot(k1x, X) < 0.0) {                       // hacia el pasado se aleja: viene del cielo
            escaped = true;
            skyDir = fromKS(normalize(-k1x));
            break;
        }

        float h = -stepSize(r, rh);
        vec3 k2x, k2p, k3x, k3p, k4x, k4p; float k2t, k3t, k4t;
        geodesic(X + 0.5 * h * k1x, P + 0.5 * h * k1p, E, a, k2x, k2p, k2t);
        geodesic(X + 0.5 * h * k2x, P + 0.5 * h * k2p, E, a, k3x, k3p, k3t);
        geodesic(X + h * k3x, P + h * k3p, E, a, k4x, k4p, k4t);
        vec3 Xold = X;
        float tOld = tAcc;
        X += (h / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
        P += (h / 6.0) * (k1p + 2.0 * k2p + 2.0 * k3p + k4p);
        tAcc += (h / 6.0) * (k1t + 2.0 * k2t + 2.0 * k3t + k4t);

        // Cruce del plano del disco (Z = 0 en Kerr-Schild, y = 0 en el shader)
        if (Xold.z * X.z < 0.0) {
            float wz = Xold.z / (Xold.z - X.z);
            vec3 Xh = mix(Xold, X, wz);
            float rd = sqrt(max(dot(Xh.xy, Xh.xy) - a * a, 0.0));
            if (rd > rin && rd < R_DISK_FADE) {
                float omega = diskOmega(rd, a);
                float g = 1.0 / (diskUt(rd, a) * (E - omega * pphi));
                float F = diskFlux(rd, a, rin) / fluxMax;
                if (turbulence) {
                    // ángulo azimutal de Kerr-Schild salientes: x + iy = (r - ia) e^{i phi}
                    float phi = atan(Xh.y * rd + a * Xh.x, Xh.x * rd - a * Xh.y);
                    F *= diskTurbulence(rd, phi, omega, tObs + mix(tOld, tAcc, wz));
                }
                float T = Tmax * pow(F, 0.25);
                float opacity = 1.0 - smoothstep(R_DISK_OUT, R_DISK_FADE, rd);
                if (!haveDebug) {
                    debug = viewMode == 1 ? divergingMap(g) : heatMap(T / Tmax);
                    raw = vec4(g, rd, 2.0, F);
                    haveDebug = true;
                }
                color += transmit * opacity * observedBlackbody(T, g, Tmax);
                transmit *= 1.0 - opacity;
                if (transmit < 0.01 || viewMode != 0) break;
            }
        }
    }

    if (viewMode == 3) {
        finalColor = haveDebug ? raw : vec4(0.0, 0.0, escaped ? 1.0 : 0.0, 0.0);
        return;
    }
    if (viewMode != 0) {
        vec3 c = haveDebug ? debug : (escaped ? vec3(0.08) : vec3(0.0));
        finalColor = vec4(c, 1.0);
        return;
    }
    if (escaped) color += transmit * renderSky(skyDir, 1.0 / E);
    finalColor = vec4(tonemap(color * exposure), 1.0);
}
