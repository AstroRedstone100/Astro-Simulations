#include "raylib.h"
#include "rlgl.h"
#include <iostream>
#include <vector>
#include <cmath>

const float G = 1.f, c = 1.f, scale = 10.f, baseh = 0.1f, minS = 0.05f, pi = 3.14159265f; const int hscreenL = 400, hscreenB = 225, screenl = 800,
screenb = 450, max_steps = 1000;
const Vector3 camerapos = {0.f, 0.f, 60.f}, up = {0.f, 1.f, 0.f};
float mass = 1.f;

struct state {
    Vector3 velocity, position;
    state(Vector3 u, Vector3 i) {
        velocity = u;
        position = i;
    }
};

float pythagore(Vector3 a) {
    return std::sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

Vector3 addVV(Vector3 a, Vector3 b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vector3 multiplyVS(Vector3 a, float m) {
    return {a.x * m, a.y * m, a.z * m};
}

float dot(Vector3 a, Vector3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

Vector3 cross(Vector3 a, Vector3 b) {
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}

state f(state ynow, float r) {
    Vector3 X = ynow.position, V = ynow.velocity, L = cross(X, V);
    Vector3 a = multiplyVS(X, dot(L, L) * -3.f * mass / (r*r*r*r*r));
    state State(a, V);
    return State;
}

state addstates(state y1, state y2) {
    state y3(addVV(y1.velocity, y2.velocity), addVV(y1.position, y2.position));
    return y3;
}

state multiplySF(state y1, float m) {
    state y2(multiplyVS(y1.velocity, m), multiplyVS(y1.position, m));
    return y2;
}

state RK4(state ynow, float h, float rnow) {
    float H = h * 0.5f;
    state k1 = f(ynow, rnow),
    
    offsetk2 = multiplySF(k1, H),
    statek2 = addstates(ynow, offsetk2);
    float rk2 = pythagore(statek2.position);
    state k2 = f(statek2, rk2),

    offsetk3 = multiplySF(k2, H),
    statek3 = addstates(ynow, offsetk3);
    float rk3 = pythagore(statek3.position);
    state k3 = f(statek3, rk3),

    offsetk4 = multiplySF(k3, h),
    statek4 = addstates(ynow, offsetk4);
    float rk4 = pythagore(statek4.position);
    state k4 = f(statek4, rk4);

    state ynext = addstates(ynow, multiplySF(addstates(addstates(k1, k4), multiplySF(addstates(k2, k3), 2.f)), h / 6.f));
    return ynext;
}

int main() {
    float r_s = 0.f;
    std::cout << "\nEnter the mass of the Black hole in 4.2 million solar masses: ";
    std::cin >> mass; r_s = 2.f * G * mass / (c*c); float R = r_s * scale;

    InitWindow(screenl, screenb, "Black Hole Simulator");
    
    Camera3D camera = {0};
    camera.position = (Vector3){0.f, 0.f, 1.f};
    camera.target = (Vector3){0.f, 0.f, 0.f};
    camera.up = up;
    camera.fovy = screenb;
    camera.projection = CAMERA_ORTHOGRAPHIC;

    Image image = LoadImage("C:/test_sfml/Background.png");
    Texture2D texture = LoadTextureFromImage(image);
    std::vector<Vector3> positions; std::vector<Color> colors;
    for (float x = 1.5f; x < screenl; x += 3.f) {
        for (float y = 1.5f; y < screenb; y += 3.f) {
            Vector3 vec = {x - hscreenL, hscreenB - y, -60.f}, vel = multiplyVS(vec, c / pythagore(vec));
            positions.push_back({x - 400.f, 225.f - y, 0.f});
            state Y(vel, camerapos); int steplim = max_steps;
            bool aboveac = false, belowac = false, change = false;
            while (true) {
                float r = pythagore(Y.position), num = r / R; num = (num > minS) ? num : minS;
                Y = RK4(Y, num * baseh, r); float rnow = pythagore(Y.position); steplim--;
                float Yguy = Y.position.y;
                if (Yguy > 1.f) {
                    aboveac = true;
                    if (belowac) {belowac = false; change = true;}
                } if (Yguy < -1.f) {
                    belowac = true;
                    if (aboveac) {aboveac = false; change = true;}
                } if (rnow / scale > 30.f) {
                    Color col = BLACK;
                    Vector3 normvel = multiplyVS(Y.velocity, 1.f/pythagore(Y.velocity));
                    float q = 0.5f + (std::atan2(normvel.z, normvel.x) / (2.f*pi)),
                    e = 0.5f - (std::asin(normvel.y) / pi);
                    int ix = (int)(q * 267) % 267, iy = (int)(e * 150) % 150;
                    ix += (ix < 0) ? 267 : 0; iy += (iy < 0) ? 150 : 0; 
                    if (ix >= 0 && ix < 267 && iy >= 0 && iy < 150) {col = GetImageColor(image, ix, iy);}
                    colors.push_back(col);
                    break;
                } else if (steplim == 0 || rnow <= R) {colors.push_back(BLACK); break;}
                else if ((rnow >= 3.f * R && rnow <= 5.f * R && std::abs(Yguy) <= 1.f) || change) {
                    colors.push_back(YELLOW);
                    break;
                }
            }
        }
    } UnloadImage(image);

    Mesh mesh = {0};
    mesh.vertexCount = positions.size();
    mesh.vertices = (float*)MemAlloc(positions.size()*3*sizeof(float));
    mesh.colors = (unsigned char*)MemAlloc(colors.size()*4*sizeof(unsigned char));
    memcpy(mesh.vertices, positions.data(), positions.size()*sizeof(Vector3));
    memcpy(mesh.colors, colors.data(), colors.size()*sizeof(Color));
    UploadMesh(&mesh, false);
    Model model = LoadModelFromMesh(mesh);

    SetTargetFPS(60);

    while (!WindowShouldClose()) {
        BeginDrawing();
            BeginMode3D(camera);
                rlEnablePointMode();
                rlSetPointSize(3.f);
                rlDisableBackfaceCulling();
                DrawModel(model, {0.f, 0.f, 0.f}, 1.f, WHITE);
                rlEnableBackfaceCulling();
                rlDisablePointMode();
            EndMode3D();
            DrawText(TextFormat("FPS: %i", GetFPS()), 2, 2, 12, WHITE); // FPS
        EndDrawing();
    }

    UnloadModel(model);
    UnloadTexture(texture);
    CloseWindow();
    return 0;
}