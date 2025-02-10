#pragma once

#include "mesh.h"
#include "colour.h"
#include "renderer.h"
#include "light.h"
#include <iostream>
#include "vec4_simd.h"
#include "matrix_simd.h"
#include <thread>
#include <mutex>
#include <vector>

// 基础2D向量类，用于屏幕空间计算
class vec2D {
public:
    float x, y;

    vec2D() { x = y = 0.f; }
    vec2D(float _x, float _y) : x(_x), y(_y) {}

#if defined(USE_SIMD_OPTIMIZATION)
    vec2D(Vec4SIMD v) {
        x = v[0];
        y = v[1];
    }
#else
    vec2D(vec4 v) {
        x = v[0];
        y = v[1];
    }
#endif

    void display() { std::cout << x << '\t' << y << std::endl; }

    vec2D operator- (vec2D& v) {
        vec2D q;
        q.x = x - v.x;
        q.y = y - v.y;
        return q;
    }
};

// 三角形渲染类，支持基础版本、SIMD优化和多线程渲染
class triangle {
private:
    // 根据编译选项选择顶点结构
#if defined(USE_SIMD_OPTIMIZATION)
    struct VertexSIMD {
        Vec4SIMD p;
        Vec4SIMD normal;
        colour rgb;
    };
    VertexSIMD v[3];
#else
    Vertex v[3];
#endif

    float area;
    colour col[3];

    // 渲染任务结构，用于多线程渲染
    struct RenderTask {
        int startY;
        int endY;
        Renderer* renderer;
        const Light* light;
        float ka;
        float kd;
    };

public:
    // 构造函数 - 支持SIMD和非SIMD版本
#if defined(USE_SIMD_OPTIMIZATION)
    triangle(const VertexSIMD& v1, const VertexSIMD& v2, const VertexSIMD& v3) {
        v[0] = v1;
        v[1] = v2;
        v[2] = v3;

        vec2D e1 = vec2D(v[1].p - v[0].p);
        vec2D e2 = vec2D(v[2].p - v[0].p);
        area = abs(e1.x * e2.y - e1.y * e2.x);
    }
#else
    triangle(const Vertex& v1, const Vertex& v2, const Vertex& v3) {
        v[0] = v1;
        v[1] = v2;
        v[2] = v3;

        vec2D e1 = vec2D(v[1].p - v[0].p);
        vec2D e2 = vec2D(v[2].p - v[0].p);
        area = abs(e1.x * e2.y - e1.y * e2.x);
    }
#endif

    // 计算重心坐标的辅助函数
    float getC(vec2D v1, vec2D v2, vec2D p) {
        vec2D e = v2 - v1;
        vec2D q = p - v1;
        return q.y * e.x - q.x * e.y;
    }

    // 计算点p的重心坐标
    bool getCoordinates(vec2D p, float& alpha, float& beta, float& gamma) {
#if defined(USE_SIMD_OPTIMIZATION)
        alpha = getC(vec2D(v[0].p), vec2D(v[1].p), p) / area;
        beta = getC(vec2D(v[1].p), vec2D(v[2].p), p) / area;
        gamma = getC(vec2D(v[2].p), vec2D(v[0].p), p) / area;
#else
        alpha = getC(vec2D(v[0].p), vec2D(v[1].p), p) / area;
        beta = getC(vec2D(v[1].p), vec2D(v[2].p), p) / area;
        gamma = getC(vec2D(v[2].p), vec2D(v[0].p), p) / area;
#endif

        if (alpha < 0.f || beta < 0.f || gamma < 0.f) return false;
        return true;
    }

    // 插值函数模板
    template <typename T>
    T interpolate(float alpha, float beta, float gamma, T a1, T a2, T a3) {
        return (a1 * alpha) + (a2 * beta) + (a3 * gamma);
    }

    // 基础版本的渲染函数
    void draw(Renderer& renderer, Light& L, float ka, float kd) {
        vec2D minV, maxV;
        getBoundsWindow(renderer.canvas, minV, maxV);

        if (area < 1.f) return;

#if defined(USE_SIMD_OPTIMIZATION)
        // SIMD优化版本的渲染循环
        for (int y = (int)(minV.y); y < (int)ceil(maxV.y); y++) {
            for (int x = (int)(minV.x); x < (int)ceil(maxV.x); x++) {
                float alpha, beta, gamma;

                if (getCoordinates(vec2D((float)x, (float)y), alpha, beta, gamma)) {
                    colour c = interpolate(beta, gamma, alpha, v[0].rgb, v[1].rgb, v[2].rgb);
                    c.clampColour();

                    float depth = interpolate(beta, gamma, alpha, v[0].p[2], v[1].p[2], v[2].p[2]);
                    Vec4SIMD normal = interpolate(beta, gamma, alpha, v[0].normal, v[1].normal, v[2].normal);
                    normal.normalise();

                    if (renderer.zbuffer(x, y) > depth && depth > 0.01f) {
                        L.omega_i.normalise();
                        float dot = max(Vec4SIMD::dot(Vec4SIMD(L.omega_i), normal), 0.0f);
                        colour a = (c * kd) * (L.L * dot + (L.ambient * ka));

                        unsigned char r, g, b;
                        a.toRGB(r, g, b);
                        renderer.canvas.draw(x, y, r, g, b);
                        renderer.zbuffer(x, y) = depth;
                        renderer.recordPixelProcessed(); // 添加像素统计
                    }
                }
            }
        }
#else
        // 基础版本的渲染循环
        for (int y = (int)(minV.y); y < (int)ceil(maxV.y); y++) {
            for (int x = (int)(minV.x); x < (int)ceil(maxV.x); x++) {
                float alpha, beta, gamma;

                if (getCoordinates(vec2D((float)x, (float)y), alpha, beta, gamma)) {
                    colour c = interpolate(beta, gamma, alpha, v[0].rgb, v[1].rgb, v[2].rgb);
                    c.clampColour();

                    float depth = interpolate(beta, gamma, alpha, v[0].p[2], v[1].p[2], v[2].p[2]);
                    vec4 normal = interpolate(beta, gamma, alpha, v[0].normal, v[1].normal, v[2].normal);
                    normal.normalise();

                    if (renderer.zbuffer(x, y) > depth && depth > 0.01f) {
                        L.omega_i.normalise();
                        float dot = max(vec4::dot(L.omega_i, normal), 0.0f);
                        colour a = (c * kd) * (L.L * dot + (L.ambient * ka));

                        unsigned char r, g, b;
                        a.toRGB(r, g, b);
                        renderer.canvas.draw(x, y, r, g, b);
                        renderer.zbuffer(x, y) = depth;
                        renderer.recordPixelProcessed(); // 添加像素统计
                    }
                }
            }
        }
#endif
    }

    // 多线程渲染的主函数
    void drawThreaded(Renderer& renderer, Light& L, float ka, float kd) {
        vec2D minV, maxV;
        getBoundsWindow(renderer.canvas, minV, maxV);

        if (area < 1.f || (maxV.y - minV.y) < 50) {
            draw(renderer, L, ka, kd);
            return;
        }

        int totalHeight = static_cast<int>(ceil(maxV.y) - floor(minV.y));
        int heightPerThread = totalHeight / MAX_THREADS;
        std::vector<std::thread> threads;

        for (int i = 0; i < MAX_THREADS; ++i) {
            RenderTask task;
            task.startY = static_cast<int>(floor(minV.y)) + i * heightPerThread;
            task.endY = (i == MAX_THREADS - 1) ?
                static_cast<int>(ceil(maxV.y)) :
                task.startY + heightPerThread;
            task.renderer = &renderer;
            task.light = &L;
            task.ka = ka;
            task.kd = kd;

            threads.emplace_back(&triangle::renderSegment, this,
                task, static_cast<int>(floor(minV.x)),
                static_cast<int>(ceil(maxV.x)));
        }

        for (auto& thread : threads) {
            thread.join();
        }
    }

    // 获取三角形在屏幕空间的边界
    void getBounds(vec2D& minV, vec2D& maxV) {
        minV = vec2D(v[0].p);
        maxV = vec2D(v[0].p);
        for (unsigned int i = 1; i < 3; i++) {
            minV.x = min(minV.x, v[i].p[0]);
            minV.y = min(minV.y, v[i].p[1]);
            maxV.x = max(maxV.x, v[i].p[0]);
            maxV.y = max(maxV.y, v[i].p[1]);
        }
    }

    // 获取经过canvas裁剪的三角形边界
    void getBoundsWindow(GamesEngineeringBase::Window& canvas, vec2D& minV, vec2D& maxV) {
        getBounds(minV, maxV);
        minV.x = max(minV.x, 0);
        minV.y = max(minV.y, 0);
        maxV.x = min(maxV.x, canvas.getWidth());
        maxV.y = min(maxV.y, canvas.getHeight());
    }

    // 调试函数：绘制三角形边界
    void drawBounds(GamesEngineeringBase::Window& canvas) {
        vec2D minV, maxV;
        getBounds(minV, maxV);

        for (int y = (int)minV.y; y < (int)maxV.y; y++) {
            for (int x = (int)minV.x; x < (int)maxV.x; x++) {
                canvas.draw(x, y, 255, 0, 0);
            }
        }
    }

    // 调试函数：显示三角形顶点信息
    void display() {
        for (unsigned int i = 0; i < 3; i++) {
            v[i].p.display();
        }
        std::cout << std::endl;
    }

private:
    // 多线程渲染中单个线程的渲染函数
    void renderSegment(RenderTask task, int minX, int maxX) {
        for (int y = task.startY; y < task.endY; ++y) {
            for (int x = minX; x < maxX; ++x) {
                float alpha, beta, gamma;

                if (getCoordinates(vec2D((float)x, (float)y), alpha, beta, gamma)) {
                    processPixel(task, x, y, alpha, beta, gamma);
                }
            }
        }
    }

    // 处理单个像素的渲染
    void processPixel(const RenderTask& task, int x, int y, float alpha, float beta, float gamma) {
#if defined(USE_SIMD_OPTIMIZATION)
        colour c = interpolate(beta, gamma, alpha, v[0].rgb, v[1].rgb, v[2].rgb);
        c.clampColour();

        float depth = interpolate(beta, gamma, alpha, v[0].p[2], v[1].p[2], v[2].p[2]);
        Vec4SIMD normal = interpolate(beta, gamma, alpha, v[0].normal, v[1].normal, v[2].normal);
        normal.normalise();

        if (depth > 0.01f) {
            
            task.light->omega_i.normalise();
            float dot = max(Vec4SIMD::dot(Vec4SIMD(task.light->omega_i), normal), 0.0f);
            colour a = (c * task.kd) * (task.light->L * dot + (task.light->ambient * task.ka));

            unsigned char r, g, b;
            a.toRGB(r, g, b);
            task.renderer->drawPixelThreadSafe(x, y, r, g, b, depth);
            task.renderer->recordPixelProcessed(); // 使用 task.renderer 而不是 renderer
        }
#else
        colour c = interpolate(beta, gamma, alpha, v[0].rgb, v[1].rgb, v[2].rgb);
        c.clampColour();

        float depth = interpolate(beta, gamma, alpha, v[0].p[2], v[1].p[2], v[2].p[2]);
        vec4 normal = interpolate(beta, gamma, alpha, v[0].normal, v[1].normal, v[2].normal);
        normal.normalise();

        if (depth > 0.01f) {
            
            vec4 light_dir = task.light->omega_i;
            light_dir.normalise();

            float dot = max(vec4::dot(light_dir, normal), 0.0f);

            // 修正颜色计算
            colour diffuse = task.light->L;
            diffuse = diffuse * dot;  // 使用已定义的 colour * float 运算符

            colour ambient = task.light->ambient;
            ambient = ambient * task.ka;  // 使用已定义的 colour * float 运算符

            colour final_color = c;
            final_color = final_color * task.kd;  // 使用已定义的 colour * float 运算符

            // 组合光照结果
            colour light_result = diffuse + ambient;  // 使用已定义的 colour + colour 运算符
            final_color = final_color * light_result;  // 使用已定义的 colour * colour 运算符

            unsigned char r, g, b;
            final_color.toRGB(r, g, b);
            task.renderer->drawPixelThreadSafe(x, y, r, g, b, depth);
            task.renderer->recordPixelProcessed(); // 使用 task.renderer 而不是 renderer
        }
#endif
    }
};