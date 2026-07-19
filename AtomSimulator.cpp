#include "raylib.h"
#include "rlgl.h"
#include <cmath>
#include <vector>
#include <random>
#include <iostream>

float render = 600.f, R = 500.f, sensitivity = 0.1f, speed = 15.f, scale = 270.f;
static std::mt19937 gen(std::random_device{}());
static std::uniform_int_distribution<> dist(-R, R);
static std::uniform_real_distribution<float> distr(0.0f, 1.0f);

Vector3 addVV(Vector3 a, Vector3 b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vector3 multiplyVS(Vector3 a, float m) {
    return {a.x * m, a.y * m, a.z * m};
}

Vector3 cross(Vector3 a, Vector3 b) {
    return {a.y*b.z - a.z*b.y, -a.x*b.z + a.z*b.x, a.x*b.y - a.y*b.x};
}

float pythagore(Vector3 a) {
    return std::sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

const std::vector<int> capacities = {0, 2, 6, 10}, subshells = {11, 21, 22, 31, 32, 41, 33}; // Number-Orbital
std::vector<int> Electronic_Configurator(int z) {
    std::vector<int> econfig; int zcopy = z;
    for (auto const& subshell : subshells) {
        int cap = capacities[subshell % 10], number = cap;
        z -= cap;
        if (z < 1) {
            number += z;
            econfig.push_back(subshell * 100 + number);
            break;
        } econfig.push_back(subshell * 100 + number);
    } if (zcopy == 24 || zcopy == 29) {
        econfig[6]++;
        econfig[5]--;
    } return econfig;
}

static const std::vector<float> maxProbs = {
    0.781212f,    // 1s
    3.112452f,    // 2s
    0.538059f,    // 2px
    0.537921f,    // 2py
    0.538032f,    // 2pz
    566.826765f,  // 3s
    17.224104f,   // 3px
    17.225414f,   // 3py
    17.224778f,   // 3pz
    447.744475f,  // 4s
    3.415680f,    // 3dxy
    0.854113f,    // 3dyz
    13.621787f,   // 3dz²
    0.853767f,    // 3dxz
    3.415693f     // 3dx²-y²
};

float normalizedchance(float r, float z, float y, float x, float a, int idx) {
    float theta = std::acos(z / r), phi = std::atan2(y, x),
    sintheta = std::sin(theta), sinphi = std::sin(phi), cosphi = std::cos(phi),
    sin2theta = sintheta*sintheta, sin2phi = sinphi*sinphi, cos2theta = z*z / (r*r), cos2phi = cosphi*cosphi, a2 = a*a,
    f = (6-a), prob = 0.f;
    if (idx == 0) {prob =  std::exp(-2.f*a);} // 1s
    else if (idx == 1) { // 2s
        float num = (2-a);
        prob = num*num*std::exp(-a);
    } else if (idx == 2) {prob = a2 * std::exp(-a) * sin2theta * cos2phi;} // 2px
    else if (idx == 3) {prob = a2 * std::exp(-a) * sin2theta * sin2phi;} // 2py
    else if (idx == 4) {prob = a2 * std::exp(-a) * cos2theta;} // 2pz
    else if (idx == 5) { // 3s
        float num = (27 - 18*a + 2*a2);
        prob = num*num * std::exp(-2.f*a/3.f);
    } else if (idx == 6) {prob = a2 * f*f * std::exp(-2.f*a/3.f) * sin2theta * cos2phi;} // 3px
    else if (idx == 7) {prob = a2 * f*f * std::exp(-2.f*a/3.f) * sin2theta * sin2phi;} // 3py
    else if (idx == 8) {prob = a2 * f*f*  std::exp(-2.f*a/3.f) * cos2theta;} // 3pz
    else if (idx == 9) { // 4s
        float num = 24.f - 18.f*a + 3.f*a2 - a2*a/6.f;
        prob = num*num * std::exp(-a/2.f);
    } else if (idx == 10) { // 3dxy
        float num = std::sin(2*phi);
        prob = a2*a2 * std::exp(-2.f*a/3.f) * sin2theta*sin2theta * num*num;
    } else if (idx == 11) {prob = a2*a2 * std::exp(-2.f*a/3.f) * sin2theta * cos2theta * sin2phi;} // 3dyz
    else if (idx == 12) { // 3dz^2
        float num = 3.f * cos2theta - 1.f;
        prob = a2*a2 * std::exp(-2.f*a/3.f) * num*num;
    } else if (idx == 13) {prob = a2*a2 * std::exp(-2.f*a / 3.f) * sin2theta * cos2theta * cos2phi;} // 3dxz
    else { // 3dx^2 - y^2
        float num = std::cos(2*phi);
        prob = a2*a2 * std::exp(-2.f*a/3.f) * sin2theta * sin2theta * num*num;
    } return prob / maxProbs[idx];
}

const std::vector<char> orbitals = {'N', 's', 'p', 'd'};
std::vector<std::vector<int>> nns = {{00}, {10}, {20, 30, 40}, {50}, {60, 70, 80}, {90}, {100, 110, 120, 130, 140}};
const std::vector<Color> colours = {WHITE, (Color){255, 255, 120, 255}, GREEN, ORANGE, (Color){0, 255, 255, 255}, RED, PURPLE};
const std::vector<int> neutrons = {0, 2, 4, 5, 6, 6, 7, 8, 10, 10, 12, 12, 14, 14, 16, 16, 18, 22, 20, 20, 24, 26, 28, 28, 30, 30,
                                   32, 30, 34, 34};

int main() {
    int k = 1;
    std::cout << "Note: This program works flawlessly only for the first 30 elements";
    std::cout << "\nEnter the atomic number Z: "; std::cin >> k; int n = neutrons[k - 1] + k;
    std::vector<int> econfig = Electronic_Configurator(k), Ints(k);
    std::cout << "\nThe electronic configuration of this element is: \n";
    
    const int pointcnt = 170000;
    std::vector<Vector3> positions; std::vector<Color> colors;
    bool flag = false; int cnt = 0;
    for (int Idx = 0; Idx < econfig.size(); Idx++) {
        int thing = econfig[Idx], subshell = thing / 100, eNum = thing % 100;
        std::cout << (subshell / 10) << orbitals[subshell % 10] << eNum << " ";
        for (int idx = 0; idx < nns[Idx].size(); idx++) {
            nns[Idx][idx] += 2; eNum -= 2;
            while (k < 0) {
                eNum++;
                nns[Idx][idx]--;
                flag = true;
            } int NUM = nns[Idx][idx]; cnt += (int)((NUM % 10) * pointcnt/ k);
            while (positions.size() < cnt) {
                float x = (float)(dist(gen)), y = (float)(dist(gen)), z = (float)(dist(gen)); Vector3 vec = {x, y, z};
                float r = pythagore(vec);
                if (r > R || r < 33.33f) {continue;}
                float chance = normalizedchance(r, z, y, x, r/scale, NUM / 10);
                if (distr(gen) > chance) {continue;}
                positions.push_back(vec);
                colors.push_back(colours[Idx]);
            } if (flag) {flag = false; break;}
        }
    } std::cout << std::endl;
    
    if (false) { // <-- Histogram here
        std::vector<int> bins(50);
        for (auto const& position : positions) {
            float radius = pythagore(position);
            bins[(int)(radius / 10.f)]++;
        } for (int idx = 0; idx < 50; idx++) {
            std::cout << idx << " --> " << bins[idx] << std::endl;
        } std::cout << std::endl;
    }

    InitWindow(800, 450, "Atom Simulator");

    Mesh mesh = {0};
    mesh.vertexCount = pointcnt;
    mesh.triangleCount = pointcnt / 3;
    mesh.vertices = (float*)MemAlloc(pointcnt*3*sizeof(float));
    mesh.colors = (unsigned char*)MemAlloc(pointcnt*4*sizeof(unsigned char));
    memcpy(mesh.vertices, positions.data(), pointcnt*sizeof(Vector3));
    memcpy(mesh.colors, colors.data(), pointcnt*sizeof(Color));
    UploadMesh(&mesh, false);
    Model model = LoadModelFromMesh(mesh);

    HideCursor();
    SetTargetFPS(60);

    Camera3D camera;
    camera.position = (Vector3){0.f, 0.f, render};
    camera.target = (Vector3){0.f, 0.f, 0.f};
    camera.up = (Vector3){0.f, 1.f, 0.f};
    camera.fovy = 90.f;
    camera.projection = CAMERA_PERSPECTIVE;

    float theta = 90.f, phi = 90.f;
    while (!WindowShouldClose()) {
        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            Vector2 mousedel = GetMouseDelta();
            theta += mousedel.y * sensitivity;
            phi -= mousedel.x * sensitivity;
            if (theta > 360.f) {theta = 0.f;} else if (theta < 0.f) {theta = 360.f;}
            if (phi > 179.f) {phi = 179.f;} else if (phi < 1.f) {phi = 1.f;}
            float thetaR = theta * DEG2RAD, phiR = phi * DEG2RAD, sinr = std::sin(thetaR);
            camera.target = addVV(camera.position, multiplyVS({sinr*std::cos(phiR), std::cos(thetaR), -sinr*std::sin(phiR)}, render));
        } if (IsKeyDown(KEY_W) || IsKeyDown(KEY_S)) {
            float sign = (IsKeyDown(KEY_W)) ? 1.f : -1.f;
            Vector3 vec = multiplyVS(addVV(camera.target, multiplyVS(camera.position, -1.f)), sign * speed / render);
            camera.position = addVV(camera.position, vec);
            camera.target = addVV(camera.target, vec);
        } else if (IsKeyDown(KEY_A) || IsKeyDown(KEY_D)) {
            float sign = (IsKeyDown(KEY_D)) ? -1.f : 1.f;
            Vector3 vec = multiplyVS(cross({0.f, 1.f, 0.f}, addVV(camera.target, multiplyVS(camera.position, -1.f))), sign * speed / render);
            camera.position = addVV(camera.position, vec);
            camera.target = addVV(camera.target, vec);
        } BeginDrawing();
            ClearBackground(BLACK);
            BeginMode3D(camera);
                DrawSphere({0.f, 0.f, 0.f}, 33.33f, BLUE); // Nucleus
                rlEnablePointMode();
                rlSetPointSize(2.f);
                rlDisableBackfaceCulling();
                DrawModel(model, {0.f, 0.f, 0.f}, 1.f, WHITE); // Electrons' superposition
                rlEnableBackfaceCulling();
                rlDisablePointMode();
            EndMode3D();
            DrawText(TextFormat("FPS: %i", GetFPS()), 2, 2, 20, WHITE); // FPS
            DrawText("+", 385, 210, 30, WHITE); // Crosshair
        EndDrawing();
    }

    CloseWindow();
    UnloadModel(model);
    return 0;
}