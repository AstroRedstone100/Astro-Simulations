#include "glad/glad.h"
#include "raylib.h"
#include "rlgl.h"
#include <vector>
#include <random>
#include <iostream>

const float side = 40.f, gap = 4.f, epsilon = 1e-3, Gap = gap*gap*2, subH = 480.f - side, subW = 800.f - side,
            side20 = 2.f*side, side15 = 1.5f*side;
const Vector2 zero = {0.f, 0.f};
const int width = 20, height = 12, fluidW = width - 2, fluidH = height - 2, vnum = 2*fluidW*fluidH - fluidH - fluidW,
          FFVnum = fluidH*fluidW - fluidH, FHnum = fluidW-1, FVnum = fluidH-1, constant11 = width+1, cnst48 = fluidH*fluidW-1,
          constant48 = (cnst48/fluidW)*2 + constant11 + cnst48, fluidSize = fluidH*fluidW;

float rho = 1.f; // Density

const Color red = {180, 4, 38, 255}, blue = {59, 76, 192, 255},
        dred = {(unsigned char)(red.r - DARKGRAY.r), (unsigned char)(red.g - DARKGRAY.g), (unsigned char)(red.b - DARKGRAY.b, 255)}, 
        dblue = {(unsigned char)(blue.r - DARKGRAY.r), (unsigned char)(blue.g - DARKGRAY.g), (unsigned char)(blue.b - DARKGRAY.b, 255)};

struct cell {
    Vector2 position;
    float pressure, divergence;
    bool block;
    std::vector<Vector2> cell_vels; // L, T, R, B
    Color color;
    cell(Vector2 pos, bool Bool) {
        position = pos;
        pressure = 0.f; divergence = 0.f;
        block = Bool;
        cell_vels = {zero, zero, zero, zero};
        color = DARKGRAY;
    }
};

float sgn(float num) {
    if (num > 0.f) {return 1.f;}
    else if (num < 0.f) {return -1.f;}
    return 0.f;
}

void DrawArrow(int x1, int y1, int x2, int y2, Color color) {
    DrawLine(x1, y1, x2, y2, color); // Vector line
    int Dx = x2 - x1, Dy = y2 - y1;
    Vector2 vertex2 = (Dx == 0.f)? (Vector2){x2 - gap, y2 - sgn(Dy)*gap} : (Vector2){x2 - sgn(Dx)*gap, y2 - gap},
    vertex3 = (Dx == 0.f)? (Vector2){x2 + gap, y2 - sgn(Dy)*gap} : (Vector2){x2 - sgn(Dx)*gap, y2 + gap};
    if (Dx < 0.f || Dy > 0.f) {DrawTriangle({(float)x2, (float)y2}, vertex3, vertex2, color);}
    else {DrawTriangle({(float)x2, (float)y2}, vertex2, vertex3, color);} // Vector Top
    DrawCircle(x1, y1, gap, color); // Vector Origin
}

void assignVelocities(std::vector<Vector2>& velocities, std::vector<cell>& cells) {
    int indx = -1;
    for (int idx = 0; idx < cells.size(); idx++) {
        if (cells[idx].block) {continue;}
        indx++;
        int xcomp = indx % fluidW, ycomp = indx / fluidW,
        r_in = xcomp + ycomp * FHnum, l_in = r_in - 1,
        b_in = FFVnum + xcomp + ycomp * fluidW, t_in = b_in - fluidW;

        if (xcomp != 0) {cells[idx].cell_vels[0] = velocities[l_in];}
        if (ycomp != 0) {cells[idx].cell_vels[1] = velocities[t_in];}
        if (xcomp != FHnum) {cells[idx].cell_vels[2] = velocities[r_in];}
        if (ycomp != FVnum) {cells[idx].cell_vels[3] = velocities[b_in];}
    }
}

float calculateDivergence(const std::vector<Vector2>& cell_vels) {
    return (cell_vels[2].x - cell_vels[0].x + cell_vels[1].y - cell_vels[3].y) / side;
}

const int range = (int)(side/2.f)-1; const float lowerBound = side/4.f,
    maxD = calculateDivergence({{-range, 0.f}, {0.f, range}, {range, 0.f}, {0.f, -range}});

void colorCells(std::vector<cell>& cells) {
    for (auto& Cell : cells) {
        float k = Cell.divergence / maxD;
        if (k > 0.f) {
            Cell.color = (Color){(unsigned char)(DARKGRAY.r + k * dred.r),
                                 (unsigned char)(DARKGRAY.g + k * dred.g), 
                                 (unsigned char)(DARKGRAY.b + k * dred.b), 255};
        } else {
            float absk = std::abs(k);
            Cell.color = (Color){(unsigned char)(DARKGRAY.r + absk * dblue.r), 
                                 (unsigned char)(DARKGRAY.g + absk * dblue.g),
                                 (unsigned char)(DARKGRAY.b + absk * dblue.b), 255};
        }
    }
}

Vector2 solveCellsFromVel(int& index) {
    if (index < FFVnum) {
        int k2 = (index/FHnum)*width + constant11 + (index%FHnum);
        return {(float)k2, (float)(k2 + 1)};
    } int k = index - FFVnum, k2 = (k/fluidW)*width + constant11 + (k%fluidW);
    return {(float)k2, (float)(k2 + width)};
    
}

float maxDivergence(const std::vector<cell>& cells) {
    float max = 0.f;
    for (auto const& Cell : cells) {
        float absDivergence = std::abs(Cell.divergence);
        if (absDivergence > max) {max = absDivergence;}
    } return max;
}

float max(const std::vector<float>& values) {
    float maximum = 0.f;
    for (auto const& value : values) {
        float absVal = std::abs(value);
        if (absVal > std::abs(maximum)) {maximum = value;}
    } return maximum;
}

void updateDivergences(std::vector<cell>& cells) {
    for (auto& Cell : cells) {
        Cell.divergence = calculateDivergence(Cell.cell_vels);
    } colorCells(cells);
}

inline Vector2 operator+(const Vector2& a, const Vector2& b) {
    return {a.x + b.x, a.y + b.y};
}

inline Vector2 operator-(const Vector2& a, const Vector2& b) {
    return {a.x - b.x, a.y - b.y};
}

inline std::vector<float> operator-(const std::vector<float>& a, const std::vector<float>& b) {
    std::vector<float> c;
    for (int idx = 0; idx < a.size(); idx++) {
        c.push_back(a[idx] - b[idx]);
    } return c;
}

const char* glslcode = R"(#version 430 core

uniform float rho;
uniform float dt;
uniform float side;
uniform int fluidwidth;
uniform int fluidheight;

layout(std430, binding = 0) buffer pressure_buffer {float pressures[];};
layout(std430, binding = 1) buffer divergence_buffer {float divergences[];};
layout(std430, binding = 2) buffer pressure_buffer2 {float new_pressures[];};
layout(local_size_x = 8, local_size_y = 8) in;

void main() {
    const int x_1 = fluidwidth - 1, y_1 = fluidheight - 1;
    uvec3 id = gl_GlobalInvocationID;
    if (id.x >= fluidwidth || id.y >= fluidheight) {return;}

    int idx = int(id.y) * fluidwidth + int(id.x);
    float sum = 0.f, self_pressure = pressures[idx];
    sum += (id.x == 0)? self_pressure : pressures[idx-1];
    sum += (id.x == x_1)? self_pressure : pressures[idx+1];
    sum += (id.y == 0)? self_pressure : pressures[idx-fluidwidth];
    sum += (id.y == y_1)? self_pressure : pressures[idx+fluidwidth];

    new_pressures[idx] = (sum - rho*side*side*divergences[idx]/dt)/4.f;
})"; // GLSL code upside

int main() {
    InitWindow(800, 480, "Fluid Simulator v1.1.0"); // Initializing all graphics
    SetTargetFPS(60); // Set FPS

    GLuint shaderObject = glCreateShader(GL_COMPUTE_SHADER); // Pass the code to the GPU
    glShaderSource(shaderObject, 1, &glslcode, NULL);
    glCompileShader(shaderObject);

    GLuint solver = glCreateProgram(); // Get the program ID
    glAttachShader(solver, shaderObject);
    glLinkProgram(solver);

    glDeleteShader(shaderObject); // Delete unnecessary shader

    static std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<> dist(-range, range); // Random initialization

    std::vector<cell> cells;
    for (float y = 0.f; y < 480.f; y += side) { // Creating cells
        for (float x = 0.f; x < 800.f; x += side) {
            bool boolean = (y == 0.f || y == subH || x == 0.f || x == subW);
            cells.push_back(cell({x, y}, boolean));
        }
    }

    std::vector<Vector2> velocities;
    while (velocities.size() < vnum) { // Assigning random velocities
        float val = 0.f;
        while (std::abs(val) < lowerBound) {val = (float)dist(gen);}
        if (velocities.size() < FFVnum) {velocities.push_back({val, 0.f});}
        else {velocities.push_back({0.f, val});}
    }

    assignVelocities(velocities, cells);
    for (auto& Cell : cells) {
        Cell.divergence = calculateDivergence(Cell.cell_vels);
    } colorCells(cells);

    float dt = 1.f;
    int selectedin = -1, i_lim = 501;
    while (!WindowShouldClose()) {
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) { // Lengthening velocities
            Vector2 mousepos = GetMousePosition();
            for (int idx = 0; idx < vnum; idx++) {
                Vector2 origin = {0.f, 0.f};
                if (idx < FFVnum) {origin = {(idx%FHnum)*side + side20, (idx/FHnum)*side + side15};}
                else {int Idx = idx - FFVnum; origin = {(Idx%fluidW)*side + side15, (Idx/fluidW)*side + side20};}
                Vector2 d = mousepos - (origin + velocities[idx]);
                if (d.x*d.x + d.y*d.y <= Gap) {selectedin = idx; break;}
            }
        } if (IsMouseButtonDown(MOUSE_BUTTON_LEFT) && selectedin != -1) {
            if (selectedin < FFVnum) {
                velocities[selectedin].x += GetMouseDelta().x;
            } else {velocities[selectedin].y += GetMouseDelta().y;}
            assignVelocities(velocities, cells);
            updateDivergences(cells);
        } else {selectedin = -1;}

        if (IsKeyPressed(KEY_D) && maxDivergence(cells) > epsilon) {
            std::cout << "\n" << maxDivergence(cells);
            std::vector<float> pressures, divergences; // Pressures and divergence initialization
            for (int index = 0; index < cells.size(); index++) {
                if (index > constant48) {break;}
                if (cells[index].block) {continue;}
                pressures.push_back(cells[index].pressure);
                divergences.push_back(cells[index].divergence);
            }

            GLuint pressure_buffer; // Exporting pressures
            glGenBuffers(1, &pressure_buffer);
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, pressure_buffer);
            glBufferData(GL_SHADER_STORAGE_BUFFER, pressures.size()*sizeof(float), pressures.data(), GL_DYNAMIC_DRAW);
            GLuint divergence_buffer; // Exporting divergences
            glGenBuffers(1, &divergence_buffer);
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, divergence_buffer);
            glBufferData(GL_SHADER_STORAGE_BUFFER, divergences.size()*sizeof(float), divergences.data(), GL_DYNAMIC_DRAW);
            std::vector<float> new_pressures(fluidSize); // New pressures vector
            GLuint pressure_buffer2;
            glGenBuffers(1, &pressure_buffer2);
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, pressure_buffer2);
            glBufferData(GL_SHADER_STORAGE_BUFFER, fluidSize*sizeof(float), new_pressures.data(), GL_DYNAMIC_DRAW);

            glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0); // Unbinding

            for (int iterations = 1; iterations < i_lim; iterations++) { // Solver loop
                glUseProgram(solver); // Use solver

                glUniform1f(glGetUniformLocation(solver, "rho"), rho); // Exporting variables
                glUniform1f(glGetUniformLocation(solver, "dt"), dt);
                glUniform1f(glGetUniformLocation(solver, "side"), side);
                glUniform1i(glGetUniformLocation(solver, "fluidwidth"), fluidW);
                glUniform1i(glGetUniformLocation(solver, "fluidheight"), fluidH);

                glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, divergence_buffer); // Linking buffers
                int num = 2*(iterations%2);
                glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2-num, pressure_buffer);
                glBindBufferBase(GL_SHADER_STORAGE_BUFFER, num, pressure_buffer2);
                glDispatchCompute((fluidW+7)/8, (fluidH+7)/8, 1); // Run calculations
                glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT); // Wait for them to finish
            }

            int finalbuffer = ((i_lim - 1)%2 == 0) ? pressure_buffer : pressure_buffer2; // Deciding final buffer
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, finalbuffer); // Binding final buffer
            glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, pressures.size()*sizeof(float), pressures.data()); // Assiging new data
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0); // Unbinding

            GLuint buffers[] = {pressure_buffer, divergence_buffer, pressure_buffer2};
            glDeleteBuffers(3, buffers); // Delete buffers
            
            for (int idx = 0; idx < pressures.size(); idx++) { // Assigning new pressures
                cells[(idx/fluidW)*2 + idx + constant11].pressure = pressures[idx];
            }
            
            for (int idx = 0; idx < vnum; idx++) { // Velocity update
                Vector2 adjacents = solveCellsFromVel(idx);
                float num = dt*(cells[(int)adjacents.y].pressure - cells[(int)adjacents.x].pressure)/(rho*side);
                if (idx < FFVnum) {velocities[idx].x -= num;}
                else {velocities[idx].y += num;}
            } assignVelocities(velocities, cells); // Re-assigning velocities to cells

            updateDivergences(cells); // Updating divergences
            std::cout << "\n" << maxDivergence(cells);
        }
        
        BeginDrawing();
            ClearBackground(BLUE);
            for (auto const& Cell : cells) { // Drawing cells
                Vector2 pos = Cell.position;
                DrawRectangleV(pos, {side, side}, ((Cell.block) ? (Color){40, 40, 40, 255} : Cell.color));
                DrawRectangleLines((int)pos.x, (int)pos.y, side, side, BLACK);
            } for (int idx = 0; idx < vnum; idx++) { // Drawing velocities
                Vector2 u = velocities[idx];
                int xcomp = (int)u.x, ycomp = (int)u.y, x = 0, y = 0;
                if (idx < FFVnum) {x = (int)((idx%FHnum)*side + side20), y = (int)((idx/FHnum)*side + side15);}
                else {int num = idx - FFVnum; x = (int)((num%fluidW)*side + side15), y = (int)((num/fluidW)*side + side20);}
                DrawArrow(x, y, x + xcomp, y + ycomp, BLUE);
            } DrawText(TextFormat("FPS: %i", GetFPS()), 2, 2, 12, WHITE); // FPS
            DrawText("[D]: Clear Divergence", 670, 466, 12, WHITE);
        EndDrawing();
    }

    CloseWindow();
    return 0;
}