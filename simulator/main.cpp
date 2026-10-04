#include "raylib.h"
#include "raymath.h"
#include "../common/KerrPhysics.h"

#include <cmath>

int main() {
    const int screenWidth = 1000;
    const int screenHeight = 800;
    InitWindow(screenWidth, screenHeight, "Real-Time Black Hole Simulator");

    ChangeDirectory(GetApplicationDirectory());
    DisableCursor();

    Camera3D camera = { 0 };
    camera.position = (Vector3){ 0.0f, 2.0f, 20.0f };
    camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy = 90.0f;
    camera.projection = CAMERA_PERSPECTIVE;

    Shader blackHoleShader = LoadShader(0, "blackhole.fs");

    int camPosLoc = GetShaderLocation(blackHoleShader, "camPos");
    int resolutionLoc = GetShaderLocation(blackHoleShader, "resolution");
    int camDirLoc = GetShaderLocation(blackHoleShader, "camDir");
    int camRightLoc = GetShaderLocation(blackHoleShader, "camRight");
    int camUpLoc = GetShaderLocation(blackHoleShader, "camUp");
    int timeLoc = GetShaderLocation(blackHoleShader, "iTime"); // Tiempo para animación
    int spinLoc = GetShaderLocation(blackHoleShader, "spin");
    int fovLoc = GetShaderLocation(blackHoleShader, "fovY");
    int tmaxLoc = GetShaderLocation(blackHoleShader, "Tmax");
    int fluxMaxLoc = GetShaderLocation(blackHoleShader, "fluxMax");
    int exposureLoc = GetShaderLocation(blackHoleShader, "exposure");
    int viewModeLoc = GetShaderLocation(blackHoleShader, "viewMode");
    int turbulenceLoc = GetShaderLocation(blackHoleShader, "turbulence");
    int timeScaleLoc = GetShaderLocation(blackHoleShader, "timeScale");

    // Parámetros físicos (unidades G = c = M = 1)
    float spin = 0.0f;          // a/M
    float lastSpin = 0.9f;      // valor al que vuelve la tecla R
    float Tmax = 6500.0f;       // K, temperatura en el pico del flujo del disco
    float exposure = 0.6f;
    int viewMode = 0;           // 0 imagen, 1 mapa de g, 2 mapa de T
    bool turbulence = true;
    float timeScale = 10.0f;    // M de tiempo por segundo
    const char* viewNames[] = { "Imagen", "Corrimiento g (azul > 1 > rojo)", "Temperatura emitida" };

    SetTargetFPS(60);

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_R)) {
            if (spin != 0.0f) { lastSpin = spin; spin = 0.0f; } else spin = lastSpin;
        }
        if (IsKeyPressed(KEY_RIGHT)) spin = fminf(spin + 0.05f, 0.95f);
        if (IsKeyPressed(KEY_LEFT)) spin = fmaxf(spin - 0.05f, -0.95f);
        if (fabsf(spin) < 1e-3f) spin = 0.0f;
        if (IsKeyPressed(KEY_UP)) exposure *= 1.25f;
        if (IsKeyPressed(KEY_DOWN)) exposure /= 1.25f;
        if (IsKeyPressed(KEY_TWO)) Tmax = fminf(Tmax * 1.25f, 1.0e6f);
        if (IsKeyPressed(KEY_ONE)) Tmax = fmaxf(Tmax / 1.25f, 1500.0f);
        if (IsKeyPressed(KEY_V)) viewMode = (viewMode + 1) % 3;
        if (IsKeyPressed(KEY_T)) turbulence = !turbulence;
        if (IsKeyPressed(KEY_PAGE_UP)) timeScale *= 2.0f;
        if (IsKeyPressed(KEY_PAGE_DOWN)) timeScale /= 2.0f;

        UpdateCamera(&camera, CAMERA_FREE);

        // La cámara es un observador que flota fuera del horizonte: no la dejamos entrar.
        kerr::Vec3 ks{ camera.position.z, camera.position.x, camera.position.y };
        double rh = kerr::horizon(spin);
        double rCam = kerr::radius(ks, spin);
        if (rCam < rh + 0.3) {
            Vector3 offset = Vector3Subtract(camera.target, camera.position);
            float push = (float)((rh + 0.3) / fmax(rCam, 1e-3));
            camera.position = Vector3Scale(camera.position, push);
            camera.target = Vector3Add(camera.position, offset);
            rCam = rh + 0.3;
        }

        Vector3 camDir = Vector3Normalize(Vector3Subtract(camera.target, camera.position));
        Vector3 camRight = Vector3Normalize(Vector3CrossProduct(camDir, camera.up));
        Vector3 camUpLocal = Vector3CrossProduct(camRight, camDir);

        float screenRes[2] = { (float)screenWidth, (float)screenHeight };
        SetShaderValue(blackHoleShader, resolutionLoc, screenRes, SHADER_UNIFORM_VEC2);

        float cameraPos[3] = { camera.position.x, camera.position.y, camera.position.z };
        SetShaderValue(blackHoleShader, camPosLoc, cameraPos, SHADER_UNIFORM_VEC3);

        float dirArr[3] = { camDir.x, camDir.y, camDir.z };
        float rightArr[3] = { camRight.x, camRight.y, camRight.z };
        float upArr[3] = { camUpLocal.x, camUpLocal.y, camUpLocal.z };
        SetShaderValue(blackHoleShader, camDirLoc, dirArr, SHADER_UNIFORM_VEC3);
        SetShaderValue(blackHoleShader, camRightLoc, rightArr, SHADER_UNIFORM_VEC3);
        SetShaderValue(blackHoleShader, camUpLoc, upArr, SHADER_UNIFORM_VEC3);

        float time = GetTime();
        SetShaderValue(blackHoleShader, timeLoc, &time, SHADER_UNIFORM_FLOAT);

        float fovY = camera.fovy * DEG2RAD;
        static float fluxMaxSpin = -2.0f, fluxMax = 1.0f;   // solo se recalcula al cambiar el spin
        if (spin != fluxMaxSpin) { fluxMax = (float)kerr::diskFluxMax(spin); fluxMaxSpin = spin; }
        int turbInt = turbulence ? 1 : 0;
        SetShaderValue(blackHoleShader, spinLoc, &spin, SHADER_UNIFORM_FLOAT);
        SetShaderValue(blackHoleShader, fovLoc, &fovY, SHADER_UNIFORM_FLOAT);
        SetShaderValue(blackHoleShader, tmaxLoc, &Tmax, SHADER_UNIFORM_FLOAT);
        SetShaderValue(blackHoleShader, fluxMaxLoc, &fluxMax, SHADER_UNIFORM_FLOAT);
        SetShaderValue(blackHoleShader, exposureLoc, &exposure, SHADER_UNIFORM_FLOAT);
        SetShaderValue(blackHoleShader, viewModeLoc, &viewMode, SHADER_UNIFORM_INT);
        SetShaderValue(blackHoleShader, turbulenceLoc, &turbInt, SHADER_UNIFORM_INT);
        SetShaderValue(blackHoleShader, timeScaleLoc, &timeScale, SHADER_UNIFORM_FLOAT);

        BeginDrawing();
            ClearBackground(BLACK);
            BeginShaderMode(blackHoleShader);
                DrawRectangle(0, 0, screenWidth, screenHeight, WHITE);
            EndShaderMode();
            DrawFPS(10, 10);
            DrawText("WASD + Raton: moverse | R: spin on/off | Izq/Der: spin | 1/2: temperatura | Arriba/Abajo: exposicion",
                     10, 40, 14, GREEN);
            DrawText("V: modo de vista | T: turbulencia | RePag/AvPag: velocidad del tiempo | ESC: salir", 10, 58, 14, GREEN);
            DrawText(TextFormat("Spin a/M = %+.2f   horizonte %.3f M   ISCO %.3f M   camara r = %.1f M",
                                spin, rh, kerr::isco(spin), rCam), 10, 84, 18, spin != 0.0f ? YELLOW : GRAY);
            DrawText(TextFormat("T max %.0f K   exposicion %.2f   turbulencia %s   tiempo x%.1f M/s",
                                Tmax, exposure, turbulence ? "si" : "no", timeScale), 10, 106, 18, GRAY);
            DrawText(TextFormat("Vista: %s", viewNames[viewMode]), 10, 128, 18, GRAY);
        EndDrawing();
    }

    UnloadShader(blackHoleShader);
    CloseWindow();
    return 0;
}
