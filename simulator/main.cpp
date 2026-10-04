#include "raylib.h"
#include "raymath.h"

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
    int kerrLoc = GetShaderLocation(blackHoleShader, "useKerr");
    int timeLoc = GetShaderLocation(blackHoleShader, "iTime"); // Tiempo para animación

    static bool rotationEnabled = false;
    SetTargetFPS(60);

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_R)) rotationEnabled = !rotationEnabled;
        UpdateCamera(&camera, CAMERA_FREE);

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

        int rotInt = rotationEnabled ? 1 : 0;
        SetShaderValue(blackHoleShader, kerrLoc, &rotInt, SHADER_UNIFORM_INT);

        float time = GetTime();
        SetShaderValue(blackHoleShader, timeLoc, &time, SHADER_UNIFORM_FLOAT);

        BeginDrawing();
            ClearBackground(BLACK);
            BeginShaderMode(blackHoleShader);
                DrawRectangle(0, 0, screenWidth, screenHeight, WHITE);
            EndShaderMode();
            DrawFPS(10, 10);
            DrawText("WASD + Mouse: Moverse | R: Alternar Rotacion (Kerr) | ESC: Salir", 10, 40, 20, GREEN);
            DrawText(TextFormat("Modo Kerr: %s", rotationEnabled ? "ACTIVADO" : "DESACTIVADO"), 10, 70, 20, rotationEnabled ? YELLOW : GRAY);
        EndDrawing();
    }

    UnloadShader(blackHoleShader);
    CloseWindow();
    return 0;
}