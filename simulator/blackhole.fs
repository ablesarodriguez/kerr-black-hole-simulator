#version 330

uniform vec3 camPos;
uniform vec2 resolution;
uniform vec3 camDir;
uniform vec3 camRight;
uniform vec3 camUp;
uniform bool useKerr;
uniform float iTime;

out vec4 finalColor;

// ---------- Ruido / estrellas ----------
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

// Rotación en 2D
vec2 rot2D(vec2 p, float a) {
    float c = cos(a);
    float s = sin(a);
    return vec2(p.x * c - p.y * s, p.x * s + p.y * c);
}

// Cielo estrellado
vec3 renderSky(vec3 dir) {
    vec2 sph = vec2(atan(dir.z, dir.x), acos(clamp(dir.y, -1.0, 1.0)));
    vec2 grid = sph * vec2(40.0, 40.0);
    vec2 id = floor(grid);
    vec2 f = fract(grid) - 0.5;

    float h = hash(id);
    float star = 0.0;
    if (h > 0.985) {
        vec2 jitter = vec2(hash(id + 3.7), hash(id + 91.3)) - 0.5;
        float d = length(f - jitter * 0.6);
        float brightness = (h - 0.985) / 0.015;
        star = smoothstep(0.12, 0.0, d) * brightness;
        float twinkle = 0.8 + 0.2 * sin(iTime * 3.0 + hash(id) * 100.0);
        star *= twinkle;
    }

    float band = pow(max(0.0, 1.0 - abs(dir.y)), 6.0) * 0.05;
    vec3 bg = vec3(0.01, 0.012, 0.02) + band * vec3(0.3, 0.28, 0.35);

    return bg + vec3(1.0, 0.97, 0.9) * star;
}

// ---------- Física ----------
vec3 get_acceleration(vec3 p, vec3 v) {
    float r2 = dot(p, p);
    if (r2 < 1.0e-6) return vec3(0.0);

    vec3 h = cross(p, v);
    float h2 = dot(h, h);

    vec3 accel = -3.0 * h2 * p / (r2 * r2 * sqrt(r2));

    if (useKerr) {
        const float spinA = 0.6;
        vec3 spinAxis = vec3(0.0, 1.0, 0.0);
        float r3 = r2 * sqrt(r2);
        accel += 2.0 * spinA * cross(spinAxis, v) / r3;
    }
    return accel;
}

void main() {
    // Cálculo de rayos UV con Jitter Temporal 
    vec2 uv = gl_FragCoord.xy / resolution.xy;
    vec2 jitter = vec2(hash(uv + iTime), hash(uv + iTime + 1.0)) - 0.5;
    uv = uv * 2.0 - 1.0 + jitter * 0.0005;
    uv.x *= resolution.x / resolution.y;

    vec3 rayDir = normalize(uv.x * camRight + uv.y * camUp + camDir);
    vec3 pos = camPos;
    vec3 vel = rayDir;

    const float R_HORIZON = 2.0;
    const float R_ISCO    = 6.0;
    const float R_DISK_OUT = 25.0;  
    const float R_DISK_FADE = 50.0; 

    bool hit_black_hole = false;
    bool hit_disk = false;
    vec3 hit_point = vec3(0.0);
    float disk_radius_hit = 0.0;
    float intensity = 0.0;
    float dopplerFactor = 1.0;
    vec3 diskTangent = vec3(0.0);
    
    float min_r = 9999.0;

    for (int i = 0; i < 1200; i++) {
        float r = length(pos);
        if (r < min_r) min_r = r;

        if (r < R_HORIZON * 1.05) { 
            hit_black_hole = true; 
            break; 
        }
        if (r > 70.0) break;

        float step_size = clamp(r * 0.04, 0.008, 0.2);
        float disk_weight = smoothstep(20.0, 6.0, r);
        step_size = mix(step_size, step_size * 0.4, disk_weight);
        step_size = clamp(step_size, 0.002, 0.2);

        vec3 k1_v = get_acceleration(pos, vel);
        vec3 k1_p = vel;
        vec3 k2_v = get_acceleration(pos + k1_p * 0.5 * step_size, vel + k1_v * 0.5 * step_size);
        vec3 k2_p = vel + k1_v * 0.5 * step_size;
        vec3 k3_v = get_acceleration(pos + k2_p * 0.5 * step_size, vel + k2_v * 0.5 * step_size);
        vec3 k3_p = vel + k2_v * 0.5 * step_size;
        vec3 k4_v = get_acceleration(pos + k3_p * step_size, vel + k3_v * step_size);
        vec3 k4_p = vel + k3_v * step_size;

        vel = normalize(vel + (step_size / 6.0) * (k1_v + 2.0*k2_v + 2.0*k3_v + k4_v));
        pos = pos + (step_size / 6.0) * (k1_p + 2.0*k2_p + 2.0*k3_p + k4_p);

        // Intersección exacta con el disco
        if (abs(vel.y) > 1.0e-6) {
            float t_plane = -pos.y / vel.y;
            if (t_plane > 0.0 && t_plane <= step_size * 1.5) {
                vec3 p_hit = pos + vel * t_plane;
                float r_xy = length(p_hit.xz);
                
                float isco_falloff = smoothstep(R_ISCO - 1.0, R_ISCO, r_xy);
                
                if (r_xy > R_ISCO - 1.0 && r_xy < R_DISK_FADE && isco_falloff > 0.01) {
                    hit_disk = true;
                    hit_point = p_hit;
                    disk_radius_hit = r_xy;
                    intensity = isco_falloff;

                    float vphi = sqrt(1.0 / disk_radius_hit);
                    diskTangent = normalize(cross(vec3(0.0, 1.0, 0.0), p_hit));

                    vec3 losDir = normalize(-vel);
                    float vLos = dot(diskTangent, losDir) * vphi;
                    float gamma = 1.0 / sqrt(max(1.0 - vphi * vphi, 1.0e-4));

                    dopplerFactor = 1.0 / (gamma * max(1.0 - vLos, 1.0e-3));
                    break;
                }
            }
        }
    }

    // ---------- Salida Final (Color cálido sin azul) ----------
    if (hit_black_hole) {
        finalColor = vec4(0.0, 0.0, 0.0, 1.0);
    } else if (hit_disk) {
        vec2 coords = hit_point.xz;
        float radius = length(coords);

        // Turbulencia
        float omega = pow(max(radius, 0.001), -1.5); 
        vec2 uv1 = rot2D(coords * 0.18, iTime * 2.0 * omega);
        float noise1 = fbm(uv1);
        vec2 uv2 = rot2D(coords * 0.22 + vec2(1.7, 9.2), iTime * 2.8 * omega + 1.2);
        float noise2 = fbm(uv2);
        float n = mix(noise1, noise2, 0.6);
        n = 0.6 + 0.4 * n;

        // FÍSICA DEL GAS
        float H = 0.25 * pow(max(radius, R_ISCO), 1.1);
        float z_abs = abs(hit_point.y);
        float density = exp(- (z_abs * z_abs) / (2.0 * H * H)); 
        float T_power = pow(max(radius, R_ISCO) / R_ISCO, -0.75);
        T_power = clamp(T_power, 0.0, 1.0);

        // MODIFICACIÓN 1: El color "caliente" ahora es amarillo dorado en lugar de blanco puro
        vec3 low = vec3(0.1, 0.0, 0.0);   // Rojo oscuro
        vec3 med = vec3(1.0, 0.2, 0.0);   // Naranja intenso
        vec3 hot = vec3(1.0, 0.9, 0.6);   // Amarillo brillante (antes era blanco con ligero azul)
        
        vec3 gas_color = low;
        if (T_power < 0.5) {
            gas_color = mix(low, med, T_power * 2.0);
        } else {
            gas_color = mix(med, hot, (T_power - 0.5) * 2.0);
        }

        // Degradado exterior
        float outer_falloff = 1.0 - smoothstep(R_DISK_OUT, R_DISK_FADE, disk_radius_hit);
        float total_intensity = intensity * outer_falloff * n * density;

        // Brillo ambiental
        float ambient_intensity = 0.05; 
        vec3 color = gas_color * (total_intensity + ambient_intensity);

        // Efecto Doppler
        float beamPower = pow(clamp(dopplerFactor, 0.05, 5.0), 3.0);
        float final_beam = 0.1 + 0.9 * beamPower; 
        color *= final_beam;

        // MODIFICACIÓN 2: El "azulado" del Doppler ha sido eliminado. 
        // Ahora al acercarse, se vuelve amarillo intenso en lugar de azul.
        color = mix(color, color * vec3(1.2, 1.0, 0.6), clamp(dopplerFactor - 1.0, 0.0, 1.0));
        color = mix(color, color * vec3(1.1, 0.4, 0.1), clamp(1.0 - dopplerFactor, 0.0, 1.0));

        // Redshift gravitacional
        float gravShift = sqrt(max(1.0 - R_HORIZON / disk_radius_hit, 1.0e-3));
        color *= gravShift; 

        // Anillo de Fotones (Dorado)
        float dist_to_photon = abs(min_r - 3.0);
        float photonRing = exp(-pow(dist_to_photon * 6.0, 2.0)) * 1.5;
        color += vec3(1.0, 0.8, 0.2) * photonRing; // Dorado puro

        finalColor = vec4(color, 1.0);
    } else {
        finalColor = vec4(renderSky(vel), 1.0);
    }
}