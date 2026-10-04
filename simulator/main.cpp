// Real-time Kerr black hole simulator: a full-screen quad shaded by blackhole.fs.
//
// Interactive by default. With --shot it renders offscreen to PNG instead (any resolution,
// optional frame sequence), which is how the images in the README are produced:
//
//   simulator --shot out.png [--size W H] [--cam X Y Z] [--target X Y Z] [--fov DEG]
//             [--spin A] [--tmax K] [--exposure E] [--view 0|1|2] [--no-turbulence]
//             [--time S] [--timescale M_PER_S] [--hud]
//             [--frames N] [--dt S] [--orbit DEG]     (writes out_0000.png, out_0001.png, ...)
//   simulator --bench N [--size W H] ...              (prints the average frame time)
#include "raylib.h"
#include "raymath.h"
#include "../common/KerrPhysics.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

static const float MAX_SPIN = 0.95f;   // the float shader is validated up to here

struct Options {
    bool offscreen = false;        // --shot or --bench
    std::string shot;              // output PNG ("" when benchmarking)
    int width = 1000, height = 800;
    Vector3 cam = { 0.0f, 4.0f, 50.0f };
    Vector3 target = { 0.0f, 0.0f, 0.0f };
    float fovDeg = 60.0f;
    float time = 0.0f;             // seconds of wall clock at the first frame
    int frames = 1;
    float dt = 1.0f / 30.0f;       // seconds between frames of a sequence
    float orbitDeg = 0.0f;         // total camera rotation about the spin axis over the sequence
    bool hud = false;
};

// Physical parameters (units G = c = M = 1)
struct Params {
    float spin = 0.0f;          // a/M
    float Tmax = 6500.0f;       // K, temperature at the peak of the disk flux
    float exposure = 0.6f;
    int viewMode = 0;           // 0 image, 1 redshift map, 2 temperature map
    bool turbulence = true;
    float timeScale = 10.0f;    // M of coordinate time per second
};

struct Uniforms {
    int camPos, resolution, camDir, camRight, camUp, time, spin, fov, tmax, fluxMax, exposure, viewMode,
        turbulence, timeScale;
};

static void usage() {
    std::puts("simulator                      interactive\n"
              "simulator --shot out.png ...   offscreen render (see the top of simulator/main.cpp)\n"
              "simulator --bench N ...        average frame time over N frames");
}

static bool parseArgs(int argc, char** argv, Options& o, Params& p) {
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        auto has = [&](int n) { return i + n < argc; };
        auto num = [&]() { return (float)std::atof(argv[++i]); };
        if (a == "--shot" && has(1)) { o.offscreen = true; o.shot = argv[++i]; }
        else if (a == "--bench" && has(1)) { o.offscreen = true; o.frames = std::atoi(argv[++i]); }
        else if (a == "--size" && has(2)) { o.width = (int)num(); o.height = (int)num(); }
        else if (a == "--cam" && has(3)) { o.cam.x = num(); o.cam.y = num(); o.cam.z = num(); }
        else if (a == "--target" && has(3)) { o.target.x = num(); o.target.y = num(); o.target.z = num(); }
        else if (a == "--fov" && has(1)) o.fovDeg = num();
        else if (a == "--time" && has(1)) o.time = num();
        else if (a == "--frames" && has(1)) o.frames = std::atoi(argv[++i]);
        else if (a == "--dt" && has(1)) o.dt = num();
        else if (a == "--orbit" && has(1)) o.orbitDeg = num();
        else if (a == "--hud") o.hud = true;
        else if (a == "--spin" && has(1)) p.spin = Clamp(num(), -MAX_SPIN, MAX_SPIN);
        else if (a == "--tmax" && has(1)) p.Tmax = num();
        else if (a == "--exposure" && has(1)) p.exposure = num();
        else if (a == "--view" && has(1)) p.viewMode = std::atoi(argv[++i]) % 3;
        else if (a == "--timescale" && has(1)) p.timeScale = num();
        else if (a == "--no-turbulence") p.turbulence = false;
        else return false;
    }
    return o.frames > 0 && o.width > 0 && o.height > 0;
}

// Sends the camera and the physical parameters to the shader and draws the full-screen quad.
static void drawBlackHole(Shader shader, const Uniforms& u, const Camera3D& camera, const Params& p,
                          float time, int width, int height) {
    Vector3 camDir = Vector3Normalize(Vector3Subtract(camera.target, camera.position));
    Vector3 camRight = Vector3Normalize(Vector3CrossProduct(camDir, camera.up));
    Vector3 camUp = Vector3CrossProduct(camRight, camDir);

    static float fluxMaxSpin = -2.0f, fluxMax = 1.0f;   // only recomputed when the spin changes
    if (p.spin != fluxMaxSpin) { fluxMax = (float)kerr::diskFluxMax(p.spin); fluxMaxSpin = p.spin; }

    float resolution[2] = { (float)width, (float)height };
    float fovY = camera.fovy * DEG2RAD;
    int turbulence = p.turbulence ? 1 : 0;
    SetShaderValue(shader, u.resolution, resolution, SHADER_UNIFORM_VEC2);
    SetShaderValue(shader, u.camPos, &camera.position, SHADER_UNIFORM_VEC3);
    SetShaderValue(shader, u.camDir, &camDir, SHADER_UNIFORM_VEC3);
    SetShaderValue(shader, u.camRight, &camRight, SHADER_UNIFORM_VEC3);
    SetShaderValue(shader, u.camUp, &camUp, SHADER_UNIFORM_VEC3);
    SetShaderValue(shader, u.time, &time, SHADER_UNIFORM_FLOAT);
    SetShaderValue(shader, u.spin, &p.spin, SHADER_UNIFORM_FLOAT);
    SetShaderValue(shader, u.fov, &fovY, SHADER_UNIFORM_FLOAT);
    SetShaderValue(shader, u.tmax, &p.Tmax, SHADER_UNIFORM_FLOAT);
    SetShaderValue(shader, u.fluxMax, &fluxMax, SHADER_UNIFORM_FLOAT);
    SetShaderValue(shader, u.exposure, &p.exposure, SHADER_UNIFORM_FLOAT);
    SetShaderValue(shader, u.viewMode, &p.viewMode, SHADER_UNIFORM_INT);
    SetShaderValue(shader, u.turbulence, &turbulence, SHADER_UNIFORM_INT);
    SetShaderValue(shader, u.timeScale, &p.timeScale, SHADER_UNIFORM_FLOAT);

    BeginShaderMode(shader);
        DrawRectangle(0, 0, width, height, WHITE);
    EndShaderMode();
}

static double cameraRadius(const Camera3D& camera, float spin) {
    return kerr::radius(kerr::Vec3{ camera.position.z, camera.position.x, camera.position.y }, spin);
}

static void drawHud(const Camera3D& camera, const Params& p, bool showKeys) {
    static const char* viewNames[] = { "Image", "Redshift g (blue > 1 > red)", "Emitted temperature" };
    int y = 10;
    if (showKeys) {
        DrawFPS(10, y);
        DrawText("WASD + mouse: move | R: spin on/off | Left/Right: spin | 1/2: temperature | Up/Down: exposure",
                 10, y + 30, 14, GREEN);
        DrawText("V: view mode | T: turbulence | PgUp/PgDn: time speed | ESC: quit", 10, y + 48, 14, GREEN);
        y += 74;
    }
    DrawText(TextFormat("Spin a/M = %+.2f   horizon %.3f M   ISCO %.3f M   camera r = %.1f M",
                        p.spin, kerr::horizon(p.spin), kerr::isco(p.spin), cameraRadius(camera, p.spin)),
             10, y, 18, p.spin != 0.0f ? YELLOW : GRAY);
    DrawText(TextFormat("T max %.0f K   exposure %.2f   turbulence %s   time x%.1f M/s",
                        p.Tmax, p.exposure, p.turbulence ? "on" : "off", p.timeScale), 10, y + 22, 18, GRAY);
    DrawText(TextFormat("View: %s", viewNames[p.viewMode]), 10, y + 44, 18, GRAY);
}

// Offscreen rendering to PNG files, or a benchmark when no output file is given.
static int runOffscreen(Shader shader, const Uniforms& u, Camera3D camera, const Options& o, const Params& p) {
    RenderTexture2D rt = LoadRenderTexture(o.width, o.height);
    std::string base = o.shot, ext = ".png";
    size_t dot = base.find_last_of('.');
    if (dot != std::string::npos) { ext = base.substr(dot); base = base.substr(0, dot); }

    double t0 = GetTime();
    for (int f = 0; f < o.frames && !WindowShouldClose(); f++) {
        float angle = o.orbitDeg * DEG2RAD * (float)f / (float)o.frames;
        camera.position = Vector3RotateByAxisAngle(o.cam, Vector3{ 0.0f, 1.0f, 0.0f }, angle);
        BeginTextureMode(rt);
            ClearBackground(BLACK);
            drawBlackHole(shader, u, camera, p, o.time + f * o.dt, o.width, o.height);
            if (o.hud) drawHud(camera, p, false);
        EndTextureMode();
        if (o.shot.empty() && f + 1 < o.frames) continue;

        Image img = LoadImageFromTexture(rt.texture);   // also waits for the GPU to finish
        if (!o.shot.empty()) {
            ImageFlipVertical(&img);
            ImageFormat(&img, PIXELFORMAT_UNCOMPRESSED_R8G8B8);
            const char* name = o.frames > 1 ? TextFormat("%s_%04d%s", base.c_str(), f, ext.c_str()) : o.shot.c_str();
            if (!ExportImage(img, name)) { UnloadImage(img); UnloadRenderTexture(rt); return 1; }
        }
        UnloadImage(img);
    }
    if (o.shot.empty()) {
        double ms = (GetTime() - t0) * 1000.0 / o.frames;
        std::printf("%dx%d  spin %+.2f  %.2f ms/frame  (%.0f fps)\n", o.width, o.height, p.spin, ms, 1000.0 / ms);
    }
    UnloadRenderTexture(rt);
    return 0;
}

int main(int argc, char** argv) {
    Options opt;
    Params p;
    if (!parseArgs(argc, argv, opt, p)) { usage(); return 1; }

    if (opt.offscreen) { SetConfigFlags(FLAG_WINDOW_HIDDEN); SetTraceLogLevel(LOG_WARNING); }
    InitWindow(opt.offscreen ? 320 : opt.width, opt.offscreen ? 240 : opt.height, "Real-Time Black Hole Simulator");

    Camera3D camera = { 0 };
    camera.position = opt.cam;
    camera.target = opt.target;
    camera.up = Vector3{ 0.0f, 1.0f, 0.0f };
    camera.fovy = opt.fovDeg;
    camera.projection = CAMERA_PERSPECTIVE;

    // The shader is copied next to the executable by CMake
    Shader shader = LoadShader(0, TextFormat("%sblackhole.fs", GetApplicationDirectory()));
    Uniforms u;
    u.camPos = GetShaderLocation(shader, "camPos");
    u.resolution = GetShaderLocation(shader, "resolution");
    u.camDir = GetShaderLocation(shader, "camDir");
    u.camRight = GetShaderLocation(shader, "camRight");
    u.camUp = GetShaderLocation(shader, "camUp");
    u.time = GetShaderLocation(shader, "iTime");
    u.spin = GetShaderLocation(shader, "spin");
    u.fov = GetShaderLocation(shader, "fovY");
    u.tmax = GetShaderLocation(shader, "Tmax");
    u.fluxMax = GetShaderLocation(shader, "fluxMax");
    u.exposure = GetShaderLocation(shader, "exposure");
    u.viewMode = GetShaderLocation(shader, "viewMode");
    u.turbulence = GetShaderLocation(shader, "turbulence");
    u.timeScale = GetShaderLocation(shader, "timeScale");

    if (opt.offscreen) {
        int status = runOffscreen(shader, u, camera, opt, p);
        UnloadShader(shader);
        CloseWindow();
        return status;
    }

    DisableCursor();
    SetTargetFPS(60);
    float lastSpin = 0.9f;      // value the R key goes back to

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_R)) {
            if (p.spin != 0.0f) { lastSpin = p.spin; p.spin = 0.0f; } else p.spin = lastSpin;
        }
        if (IsKeyPressed(KEY_RIGHT)) p.spin = fminf(p.spin + 0.05f, MAX_SPIN);
        if (IsKeyPressed(KEY_LEFT)) p.spin = fmaxf(p.spin - 0.05f, -MAX_SPIN);
        if (fabsf(p.spin) < 1e-3f) p.spin = 0.0f;
        if (IsKeyPressed(KEY_UP)) p.exposure *= 1.25f;
        if (IsKeyPressed(KEY_DOWN)) p.exposure /= 1.25f;
        if (IsKeyPressed(KEY_TWO)) p.Tmax = fminf(p.Tmax * 1.25f, 1.0e6f);
        if (IsKeyPressed(KEY_ONE)) p.Tmax = fmaxf(p.Tmax / 1.25f, 1500.0f);
        if (IsKeyPressed(KEY_V)) p.viewMode = (p.viewMode + 1) % 3;
        if (IsKeyPressed(KEY_T)) p.turbulence = !p.turbulence;
        if (IsKeyPressed(KEY_PAGE_UP)) p.timeScale *= 2.0f;
        if (IsKeyPressed(KEY_PAGE_DOWN)) p.timeScale /= 2.0f;

        UpdateCamera(&camera, CAMERA_FREE);

        // The camera is an observer hovering outside the horizon: keep it from going in.
        double rMin = kerr::horizon(p.spin) + 0.3, rCam = cameraRadius(camera, p.spin);
        if (rCam < rMin) {
            Vector3 offset = Vector3Subtract(camera.target, camera.position);
            camera.position = Vector3Scale(camera.position, (float)(rMin / fmax(rCam, 1e-3)));
            camera.target = Vector3Add(camera.position, offset);
        }

        BeginDrawing();
            ClearBackground(BLACK);
            drawBlackHole(shader, u, camera, p, (float)GetTime() + opt.time, opt.width, opt.height);
            drawHud(camera, p, true);
        EndDrawing();
    }

    UnloadShader(shader);
    CloseWindow();
    return 0;
}
