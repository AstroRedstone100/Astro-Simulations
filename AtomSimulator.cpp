#include "raylib.h"
#include "rlgl.h"
#include <unordered_set>
#include <cmath>
#include <vector>
#include <random>
#include <iostream>

const float render = 600.f, R = 500.f, sensitivity = 0.1f, scale = 270.f, r = 6.6f, r_ = 7.2f, A = 10.f;
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

// These hardcoded values are evaluated by ChatGPT by conducting a numerical search on the formulas
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

int getuniqueid(Vector3 a) {
    return (int)(a.x * 1e9) + (int)(a.y * 1e6) + (int)(a.z * 1e3);
}

const std::vector<float> angles = {90.0f, 45.0f, 135.0f, 180.0f, 22.5f, 112.5f, 67.5f, 157.5f, 33.75f, 123.75f, 56.25f, 146.25f, 11.25f,
    101.25f, 78.75f, 168.75f, 39.375f, 129.375f, 50.625f, 140.625f, 16.875f, 106.875f, 73.125f, 163.125f, 5.625f, 95.625f, 84.375f, 174.375f};

struct nucleon {
    Color color;
    Vector3 position;
    nucleon(Color colour, Vector3 pos) {
        color = colour;
        position = pos;
    }
};

int main() {
    int k = 1; float R_ = 0.f;
    std::cout << "Note: This program works flawlessly only for the first 30 elements";
    std::cout << "\nEnter the atomic number Z: "; std::cin >> k; int M = neutrons[k - 1], n = M + k;
    std::vector<int> econfig = Electronic_Configurator(k), Ints(k);
    std::cout << "\nThe electronic configuration of this element is: \n";
    
    // My own nucleus shape generation algorithm
    std::unordered_set<int> cabinet;
    std::vector<Vector3> directions;
    float y = (n % 2 == 0) ? 0.f : r; R_ = r * std::cbrt((float)n);
    int x = std::ceil((R_ - y) / (2.f * r)), cnt0 = 0, cnt1 = 0, cnt2 = 0, numbr = (int)(y/r);
    float angle1 = 0.f, angle2 = 0.f;
    Vector3 dir;
    while (directions.size() * x + numbr < n) {
        if (cnt0 == 0) {
            float angle1R = angle1 * DEG2RAD, angle2R = angle2 * DEG2RAD, sinr = std::sin(angle1R);
            dir = {sinr*std::cos(angle2R), std::cos(angle1R), -sinr*std::sin(angle2R)};
        } else {
            dir = multiplyVS(dir, -1.f);
            cnt1++; if (cnt1 > 27) {cnt1 = 0;}
            cnt2 += (cnt1 % 2 == 1) ? 1 : 0; if (cnt2 > 27) {cnt2 = 0;}
            angle1 = angles[cnt1]; angle2 = angles[cnt2];
        } cnt0 = 1 - cnt0;
        int id = getuniqueid(dir);
        if (cabinet.count(id) == 0) {
            cabinet.insert(id);
            directions.push_back(dir);
        }
    } std::vector<nucleon> nucleons;
    if (numbr != 0) {nucleon Nucleon(BLUE, {0.f, 0.f, 0.f}); nucleons.push_back(Nucleon);}
    bool breaks = false;
    for (auto const& Dir : directions) {
        for (int m = 1; m <= x; m++) {
            Color COLOR = BLUE; if (numbr % 2 == 1 && M != 0) {COLOR = RED; M--;}
            nucleon Nucleon(COLOR, multiplyVS(Dir, r*(2.f*m - 1.f) + y));
            nucleons.push_back(Nucleon);
            numbr++;
            if (numbr == n) {breaks = true; break;}
        } if (breaks) {
            breaks = false;
            break;
        }
    }

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
                if (r > R || r < R_) {continue;}
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

    float theta = 90.f, phi = 90.f, speed = 0.f;
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
            float sign = (IsKeyDown(KEY_W)) ? 1.f : -1.f; speed += GetFrameTime() * A;
            Vector3 vec = multiplyVS(addVV(camera.target, multiplyVS(camera.position, -1.f)), sign * speed / render);
            camera.position = addVV(camera.position, vec);
            camera.target = addVV(camera.target, vec);
        } else if (IsKeyDown(KEY_A) || IsKeyDown(KEY_D)) {
            float sign = (IsKeyDown(KEY_D)) ? -1.f : 1.f; speed += GetFrameTime() * A;
            Vector3 vec = multiplyVS(cross({0.f, 1.f, 0.f}, addVV(camera.target, multiplyVS(camera.position, -1.f))), sign * speed / render);
            camera.position = addVV(camera.position, vec);
            camera.target = addVV(camera.target, vec);
        } else {speed = 0.f;}
        BeginDrawing();
            ClearBackground(BLACK);
            BeginMode3D(camera);
                for (auto const& Nucleon : nucleons) {
                    DrawSphere(Nucleon.position, r_, Nucleon.color);  // Nucleons
                } rlEnablePointMode();
                rlSetPointSize(1.1f);
                rlDisableBackfaceCulling();
                DrawModel(model, {0.f, 0.f, 0.f}, 1.f, WHITE); // Electrons' superposition
                rlEnableBackfaceCulling();
                rlDisablePointMode();
            EndMode3D();
            DrawText(TextFormat("FPS: %i", GetFPS()), 2, 2, 20, WHITE); // FPS
            DrawText("+", 385, 210, 30, WHITE); // Crosshair
            DrawText(TextFormat("Speed: %.0f px/s", speed), 2, 428, 20, WHITE); // Speed
        EndDrawing();
    }

    CloseWindow();
    UnloadModel(model);
    return 0;
}