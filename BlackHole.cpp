#include "raylib.h"
#include <cmath>
#include <vector>

Vector3 addV(Vector3 a, Vector3 b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vector3 multiplyV(Vector3 a, float m) {
    return {a.x * m, a.y * m, a.z * m};
}

float pythagore(Vector3 v) {
    return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

float dot(Vector3 a, Vector3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

Vector3 cross(Vector3 a, Vector3 b) {
    return {a.y*b.z - a.z*b.y, b.x*a.z - a.x*b.z, a.x*b.y - a.y*b.x};
}

const Vector3 camerapos = {0.f, -80.f, 900.f}, bhpos = {0.f, 0.f, 0.f};
const float G = 1.f, c = 1.f, mass = 1.f, scale = 35.f, rs = 2.f*G*mass*scale/(c*c), z = 318.f, p = 3.f * rs, q = -p/2.f, u = p*3.f;

struct state {
    Vector3 position, velocity;
    state(Vector3 pos, Vector3 vel) {
        position = pos;
        velocity = vel;
    }
};

state addS(state a, state b) {
    return state(addV(a.position, b.position), addV(a.velocity, b.velocity));
}

state multiplyS(state a, float m) {
    return state(multiplyV(a.position, m), multiplyV(a.velocity, m));
}

state f(state ynow) {
    Vector3 position = ynow.position, velocity = ynow.velocity, L = cross(position, velocity);
    float r = pythagore(position);
    return state(velocity, multiplyV(position, q * dot(L, L) / (r*r*r*r*r)));
}

state rk4(state ynow, float h) {
    float hh = h/2.f;
    state k1 = f(ynow), k2 = f(addS(ynow, multiplyS(k1, hh))),
    k3 = f(addS(ynow, multiplyS(k2, hh))), k4 = f(addS(ynow, multiplyS(k3, h)));
    return addS(ynow, multiplyS(addS(addS(k1, k4), multiplyS(addS(k2, k3), 2.f)), hh/3.f));
}

float sgn(float num) {
    if (num == 0.f) {return 0.f;}
    else if (num > 0.f) {return 1.f;}
    return -1.f;
}

int main() {
    InitWindow(800, 450, "Black Hole Simulator v2.2");
    SetTargetFPS(60);

    std::vector<Color> pixels;
    for (float y = 225.f; y > -225.f; y -= 1.f) {
        for (float x = -400.f; x < 400; x += 1.f) {
            Vector3 notdir = {x, y, z};
            state y(camerapos, multiplyV(notdir, -c / pythagore(notdir)));
            while (true) {
                Vector3 pos = y.position;
                float r = pythagore(pos), Y = pos.y;
                if (r < rs || r > 2500.f) {pixels.push_back(BLACK); break;}
                y = rk4(y, std::max(0.05f, (r / rs) - 0.95f));
                if (r >= p && r <= u && (std::abs(Y) <= 1.f || sgn(Y) != sgn(y.position.y))) {pixels.push_back({255, 120, 20, 255}); break;}
            }
        }
    }

    Image image = {0};
    image.data = pixels.data();
    image.width = 800;
    image.height = 450;
    image.mipmaps = 1;
    image.format = PIXELFORMAT_UNCOMPRESSED_R8G8B8A8;

    Texture2D texture = LoadTextureFromImage(image);
    while (!WindowShouldClose()) {
        BeginDrawing();
            ClearBackground(BLUE);
            DrawTexture(texture, 0, 0, WHITE);
            DrawText(TextFormat("FPS: %i", GetFPS()), 2, 2, 12, WHITE);
        EndDrawing();
    }

    UnloadTexture(texture);
    CloseWindow();
    return 0;
}