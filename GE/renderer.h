#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include "GamesEngineeringBase.h"
#include "zbuffer.h"
#include "matrix.h"
#include "matrix_simd.h"
#include "vec4_simd.h"

#include <sstream>  // 用于std::stringstream
#include <string>   // 用于std::string
#include <iomanip>  // 用于格式化输出

// The `Renderer` class handles rendering operations, including managing the
// Z-buffer, canvas, and perspective transformations for a 3D scene.
class Renderer {
    float fov = 90.0f * M_PI / 180.0f; // Field of view in radians (converted from degrees)
    float aspect = 4.0f / 3.0f;        // Aspect ratio of the canvas (width/height)
    float n = 0.1f;                    // Near clipping plane distance
    float f = 100.0f;                  // Far clipping plane distance
    // 性能统计数据
    struct PerformanceMetrics {
        double frameTime;
        int triangleCount;
        int pixelsProcessed;
    } metrics;
public:
    Zbuffer<float> zbuffer;                  // Z-buffer for depth management
    GamesEngineeringBase::Window canvas;     // Canvas for rendering the scene
    

#if defined(USE_SIMD_OPTIMIZATION)
    MatrixSIMD perspective;
#else
    matrix perspective;
#endif


    // Constructor initializes the canvas, Z-buffer, and perspective projection matrix.
    Renderer() {
        canvas.create(1024, 768, "Raster");  // Create a canvas with specified dimensions and title
        zbuffer.create(1024, 768);           // Initialize the Z-buffer with the same dimensions
#if defined(USE_SIMD_OPTIMIZATION)
        perspective = MatrixSIMD::makePerspective(fov, aspect, n, f);
#else
        perspective = matrix::makePerspective(fov, aspect, n, f);
#endif; // Set up the perspective matrix
        resetMetrics();
    }

    // Clears the canvas and resets the Z-buffer.
    void clear() {
        canvas.clear();  // Clear the canvas (sets all pixels to the background color)
        zbuffer.clear(); // Reset the Z-buffer to the farthest depth
        resetMetrics();
    }

    // Presents the current canvas frame to the display.
    void present() {
        canvas.present(); // Display the rendered frame
        updatePerformanceMetrics();
    }
    // 获取性能指标
    const PerformanceMetrics& getPerformanceMetrics() const {
        return metrics;
    }
    // 透视变换处理
#if defined(USE_SIMD_OPTIMIZATION)
    Vec4SIMD applyPerspective(const Vec4SIMD& vertex) {
        Vec4SIMD transformed = perspective * vertex;
        transformed.divideW();
        return transformed;
    }
#else
    vec4 applyPerspective(const vec4& vertex) {
        vec4 transformed = perspective * vertex;
        transformed.divideW();
        return transformed;
    }
#endif

    // 视口变换
    template<typename VecType>
    void applyViewportTransform(VecType& vertex) {
        vertex[0] = (vertex[0] + 1.f) * 0.5f * static_cast<float>(canvas.getWidth());
        vertex[1] = (vertex[1] + 1.f) * 0.5f * static_cast<float>(canvas.getHeight());
        vertex[1] = canvas.getHeight() - vertex[1];
    }

private:
    void resetMetrics() {
        metrics = PerformanceMetrics();
    }

    void updatePerformanceMetrics() {
        // 这里可以添加帧时间计算等性能统计
        static auto lastFrameTime = std::chrono::high_resolution_clock::now();
        auto currentTime = std::chrono::high_resolution_clock::now();
        metrics.frameTime = std::chrono::duration<double, std::milli>(currentTime - lastFrameTime).count();
        lastFrameTime = currentTime;
    }

public:
    // 性能分析功能
    void startFrameMetrics() {
        resetMetrics();
    }

    void recordTriangleProcessed() {
        metrics.triangleCount++;
    }

    void recordPixelProcessed() {
        metrics.pixelsProcessed++;
    }

    // 获取性能报告
    std::string getPerformanceReport() const {
        std::stringstream ss;
        ss << "Performance Report:\n";
        ss << "Frame Time: " << metrics.frameTime << "ms\n";
        ss << "Triangles Processed: " << metrics.triangleCount << "\n";
        ss << "Pixels Processed: " << metrics.pixelsProcessed << "\n";
        ss << "Average Pixels per Triangle: "
            << (metrics.triangleCount > 0 ? metrics.pixelsProcessed / metrics.triangleCount : 0) << "\n";
        return ss.str();
    }
};

