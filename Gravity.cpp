#include "raylib.h"
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
// 0.15 AI uses, cd "c:\test_sfml\" ; if ($?) { g++ -I C:\\test_sfml\\eigen-5.0.1 Gravity.cpp -I C:\raylib\include -L C:\raylib\lib -lraylib -lopengl32 -lgdi32 -lwinmm -o output } ; if ($?) { .\output }
float M_O = 1.989e30f, R_O = 1.16e4f, view_r = 800.f, dx = 17.5f, sensitivity = 0.01f, epsilon = 0.001f, G = 6.6743e-11f,
      avg_FT = 0.0166667f, R_sun = 6.96e5f;

struct body {
    Color color;
    float radius, mass;
    Vector3 position, v_0, acting_a, W;
    body(Color colour, float r, float m, Vector3 pos, Vector3 u, Vector3 W_0) {
        color = colour;
        radius = r / R_O; mass = m / M_O;
        position = pos; v_0 = u; W = W_0;
        acting_a = {0.f, 0.f, 0.f};
    }
};

struct point {
    Vector3 position;
    bool tobedrawn;
    point(Vector3 Pos) {
        position = Pos;
        tobedrawn = true;
    }
};

Vector3 multiplyVS(Vector3 vec, float scalar) {
    return {vec.x * scalar, vec.y * scalar, vec.z * scalar};
}

Vector3 addVV(Vector3 vec1, Vector3 vec2) {
    return {vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z};
}

float pythagore(Vector3 vec) {
    return std::sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
}

float cforce(float m1, float m2, float Radius) {
    float Force = M_O * G * m1 * m2 / std::pow(Radius * R_O * 1000.f, 2);
    return Force;
}

double dot_product(Vector3 vec1, Vector3 vec2) {
    float X = vec1.x * vec2.x, Y = vec1.y * vec2.y, Z = vec1.z * vec2.z;
    return (X + Y + Z);
}

Vector3 cross_product(Vector3 vec1, Vector3 vec2) {
    Vector3 answer = {vec1.y*vec2.z - vec1.z*vec2.y,
                      vec1.z*vec2.x - vec1.x*vec2.z,
                      vec1.x*vec2.y - vec1.y*vec2.x};
    return answer;
}

int sgn(double x86) {
    return ((x86 >= 0) ? 1 : -1);
}

bool hasgoodroots(double p, double q, double r, double s, double t) {
    if (std::abs(p) < 1e-12) return false; 
    bool ret = true;
    double b = q/p, c = r/p, d = s/p, e = t/p,
    P = 8*c - 3*b*b,
    D = 64*d - 16*b*c + 3*b*b*b,
    dis = 256*e*e*e - 192*b*d*e*e - 128*c*c*e*e + 144*c*d*d*e 
               - 27*d*d*d*d + 144*b*b*c*e*e - 6*b*b*d*d*e - 80*b*c*c*d*e 
               + 18*b*c*d*d*d + 16*c*c*c*c*e - 4*c*c*c*d*d - 27*b*b*b*b*e*e 
               + 18*b*b*b*c*d*e - 4*b*b*b*d*d*d - 4*b*b*c*c*c*e + b*b*c*c*d*d;
    if (dis > 0) {
        if (P >= 0.0 || D >= 0.0) {ret = false;}
    } if (ret && e > 0.0) {
        int cnt = 0;
        std::vector<double> coeffs = {1.0, b, c, d, e};
        double last_nonzero = 1.0;
        for (int idx = 1; idx < 5; idx++) {
            if (coeffs[idx] != 0.0) {
                if ((coeffs[idx] > 0 && last_nonzero < 0) || (coeffs[idx] < 0 && last_nonzero > 0)) {
                    cnt++;
                } last_nonzero = coeffs[idx];
            }
        } if (cnt == 0) {ret = false;}
    } return ret;
}

int main() {
    InitWindow(800, 450, "Gravity Simulator");
    SetTargetFPS(60);
    HideCursor();

    Camera3D camera;
    camera.position = (Vector3){0.f, 500.f, 400.f};
    camera.target = camera.position; camera.target.z = camera.position.z - view_r;
    camera.up = (Vector3){0.f, 1.f, 0.f};
    camera.fovy = 90.f;
    camera.projection = CAMERA_PERSPECTIVE;

    body earth(WHITE, R_sun, M_O, {0.f, 200.f, 0.f}, {0.f, 0.f, 0.f}, {5.f, 0.f, 0.f}); // NOTHING TO SCALE HERE
    body moon(YELLOW, 0.272f * R_sun, M_O / 81.f, {781.0249676f, 200.f, 0.f}, {0.f, 0.f, -50.f}, {0.f, 0.f, -5.fk}); // SAME HERE
    std::vector<body> bodies = {earth, moon};

    std::vector<point> points;
    for (float z = -2000.f; z < 2100.f; z += 100.f) {
        for (float x = -2000.f; x < 2100.f; x += 100.f) {
            point Point({x, 0.f, z});
            points.push_back(Point);
        }
    }
    float textsize = 30.f;
    Vector2 size = MeasureTextEx(GetFontDefault(), "+", textsize, 0.f),
    imp_vec = {400.f - size.x / 2.f, 225.f - size.y / 2.f};
    float theta = 90.f, phi = 90.f;
    while (!WindowShouldClose()) {
        for (auto& Point : points) {
            Point.position.y = 0.f; Point.tobedrawn = true;
            std::vector<Vector3> all;
            for (auto const& Body : bodies) {
                Vector3 forusev = Body.position; forusev.y = 0.f;
                Vector3 disVec = addVV(forusev, multiplyVS(Point.position, -1.f));
                disVec = addVV(disVec, {epsilon, epsilon, epsilon});
                float dist = pythagore(disVec);
                Vector3 massVec = multiplyVS(disVec, Body.mass / (dist*dist*dist));
                all.push_back(massVec);
            } Vector3 totalVec = {0.f, 0.f, 0.f};
            for (auto const& Vec : all) {
                totalVec = addVV(totalVec, Vec);
            } float depth = -pythagore(totalVec) * 500000.f;
            if (depth < -33000.f) {Point.tobedrawn = false;}
            else {Point.position.y = depth;}
        }
        
        std::vector<body> bodycopy = bodies;
        float dt = (float)GetFrameTime();
        std::vector<std::vector<Vector3>> forcess;
        for (int idx = 0; idx < bodies.size(); idx++) {
            std::vector<Vector3> forces;
            for (int idx2 = idx + 1; idx2 < bodies.size(); idx2++) {
                Vector3 disVec = addVV(bodies[idx2].position, multiplyVS(bodies[idx].position, -1.f));
                float R = pythagore(disVec), f_grav = cforce(bodies[idx].mass, bodies[idx2].mass, R);
                Vector3 forceV = multiplyVS(disVec, f_grav / R);
                forces.push_back(forceV);
            } forcess.push_back(forces);
            if (idx != 0) {
                for (int num = 0; num < idx; num++) {
                    forces.push_back(multiplyVS(forcess[idx - num - 1][num], -1.f));
                }
            } Vector3 net_force = {0.f, 0.f, 0.f};
            for (auto const& Force : forces) {
                net_force = addVV(net_force, Force);
            } net_force = multiplyVS(net_force, 4.8e7f);
            Vector3 acceleration = multiplyVS(net_force, 1/bodies[idx].mass),
            pos = addVV(multiplyVS(acceleration, dt*dt/2), multiplyVS(bodies[idx].v_0, dt * R_O * 1000.f)),
            velocity = addVV(multiplyVS(acceleration, dt), multiplyVS(bodies[idx].v_0, R_O * 1000.f));
            bodies[idx].acting_a = multiplyVS(acceleration, 1/(R_O * 1000.f));

            bodies[idx].v_0 = multiplyVS(velocity, 1/(R_O * 1000.f));
            bodies[idx].position = addVV(bodies[idx].position, multiplyVS(pos, 1/(R_O * 1000.f)));
        }

        while (true) {
            int colcnt = 0;
            for (int idx = 0; idx < bodies.size(); idx++) {
                for (int idx2 = idx + 1; idx2 < bodies.size(); idx2++) {
                    bool collision = false;
                    Vector3 pos1 = bodycopy[idx].position, pos2 = bodycopy[idx2].position,
                    dis_n = addVV(bodies[idx].position, multiplyVS(bodies[idx2].position, -1.f));
                    float r_n = pythagore(dis_n), r_sum = bodies[idx].radius + bodies[idx2].radius;
                    if (r_n <= r_sum/1.5f) {collision = true;} else {
                        Vector3 dis_0 = addVV(bodycopy[idx].position, multiplyVS(bodycopy[idx2].position, -1.f));
                        float r_0 = pythagore(dis_0);

                        if (r_0 >= r_n) {continue;}

                        Vector3 DA = addVV(bodies[idx].acting_a, multiplyVS(bodies[idx2].position, -1.f)),
                        DV_0 = addVV(bodycopy[idx].v_0, multiplyVS(bodycopy[idx2].v_0, -1.f)),
                        DP_0 = addVV(bodycopy[idx].position, multiplyVS(bodycopy[idx2].position, -1.f));
                        double A = 0.25 * dot_product(DA, DA), B = dot_product(DA, DV_0), C = dot_product(DV_0, DV_0) + dot_product(DP_0, DA),
                              D = 2.0 * dot_product(DV_0, DP_0), E = dot_product(DP_0, DP_0) - r_sum*r_sum/2.25;
                        collision = hasgoodroots(A, B, C, D, E);
                    }
                
                    if (collision) {
                        colcnt++;
                        float m1 = bodies[idx].mass, m2 = bodies[idx2].mass, mu = m1 * m2 / (m1 + m2),
                              I = 5.f * (m1 + m2) * r_sum*r_sum / 18.f;

                        Vector3 v1 = bodycopy[idx].v_0, v2 = bodycopy[idx2].v_0, deltaV = addVV(v1, multiplyVS(v2, -1.f)),
                                L = multiplyVS(cross_product(addVV(pos1, multiplyVS(pos2, -1.f)), deltaV), mu), w = multiplyVS(L, 1/I),
                                v_resultant = multiplyVS(addVV(multiplyVS(v1, m1), multiplyVS(v2, m2)), 1/(m1 + m2));

                        Color c1 = bodies[idx].color, c2 = bodies[idx2].color;
                        body newbody({(unsigned char)((c1.r + c2.r)/2), (unsigned char)((c1.g + c2.g)/2), (unsigned char)((c1.b + c2.b)/2), 255}, 5.f * r_sum * R_O/6.f,
                             (m1 + m2) * M_O, multiplyVS(addVV(bodycopy[idx].position, bodycopy[idx2].position), 0.5f), v_resultant, w);

                        bodies.erase(bodies.begin() + idx2);
                        bodies[idx] = newbody;
                        bodycopy.erase(bodycopy.begin() + idx2);
                        bodycopy[idx] = newbody;
                    }
                }
            } if (colcnt == 0) {break;}
        }

        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            Vector2 mousedel = GetMouseDelta();
            theta += mousedel.x * sensitivity;
            phi += mousedel.y * sensitivity;
            if (theta > 179.f) {theta = 179.f;} else if (theta < 1.f) {theta = 1.f;}
            if (phi > 179.f) {phi = 179.f;} else if (phi < 1.f) {phi = 1.f;}

            camera.target = camera.position;
            camera.target.z -= std::sin(theta) * std::sin(phi) * view_r;
            camera.target.x += std::cos(theta) * std::sin(phi) * view_r;
            camera.target.y += std::cos(phi) * view_r;

        } if (IsKeyDown(KEY_W) || IsKeyDown(KEY_S)) {
            int sign = 1; if (IsKeyDown(KEY_S)) {sign = -1;}

            camera.position.z -= sign * std::sin(theta) * std::sin(phi) * dx;
            camera.position.x += sign * std::cos(theta) * std::sin(phi) * dx;
            camera.position.y += sign * std::cos(phi) * dx;

            camera.target.z -= sign * std::sin(theta) * std::sin(phi) * dx;
            camera.target.x += sign * std::cos(theta) * std::sin(phi) * dx;
            camera.target.y += sign * std::cos(phi) * dx;
        } else if (IsKeyDown(KEY_A) || IsKeyDown(KEY_D)) {
            int sign = -1; if (IsKeyDown(KEY_D)) {sign = 1;}
            camera.position.x += sign * dx;
            camera.target.x += sign * dx;
        }
        BeginDrawing();
            ClearBackground(BLACK);
            BeginMode3D(camera);
                for (auto const& Body : bodies) {
                    DrawSphere(Body.position, Body.radius, Body.color);
                } for (int idx = 0; idx < points.size(); idx++) {
                    if (!points[idx].tobedrawn) {continue;}
                    if (idx + 41 < points.size()) {
                        if (!points[idx + 41].tobedrawn) {continue;}
                        DrawLine3D(points[idx].position,
                        points[idx + 41].position, DARKGRAY);}
                    if (idx + 1 < points.size() && (idx + 1) % 41 != 0) {
                        if (!points[idx + 1].tobedrawn) {continue;}
                        DrawLine3D(points[idx].position,
                        points[idx + 1].position, DARKGRAY);}
                }
            EndMode3D();
            DrawText((std::to_string(GetFPS()) + " FPS").c_str(), 3, 3, 34, GREEN);
            DrawTextPro(GetFontDefault(), "+", imp_vec, {0.f, 0.f}, 0.f, textsize, 0.f, WHITE); // Crosshair
        EndDrawing();
    }

    CloseWindow();
    return 0;
}