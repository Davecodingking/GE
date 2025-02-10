#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include "GamesEngineeringBase.h"
#include "zbuffer.h"
#include "matrix.h"
#include "matrix_simd.h"
#include "vec4_simd.h"
#include <thread>  // 添加线程支持
#include <vector>
#include <mutex>   // 添加互斥锁支持
#include <atomic>  // 添加原子操作支持
#include <sstream>
#include <string>
#include <iomanip>

// 定义线程数量的常量，可以根据CPU核心数调整
const int MAX_THREADS = 11;  // 根据作业要求，最多使用11个线程

class Renderer {
private:
    float fov = 90.0f * M_PI / 180.0f;
    float aspect = 4.0f / 3.0f;
    float n = 0.1f;
    float f = 100.0f;

    struct PerformanceMetrics {
        double frameTime;
        std::atomic<int> triangleCount;     // 改为原子类型,因为三角形也可能在多线程中处理
        std::atomic<int> pixelsProcessed;
        std::atomic<int> activeThreads;

        // 添加构造函数
        PerformanceMetrics() :
            frameTime(0.0),
            triangleCount(0),
            pixelsProcessed(0),
            activeThreads(0) {
        }

        // 添加自定义赋值操作
        void reset() {
            frameTime = 0.0;
            triangleCount = 0;
            pixelsProcessed = 0;
            activeThreads.store(0);
        }
    } metrics;

    // 多线程渲染需要的互斥锁
    std::mutex zbufferMutex;
    std::vector<std::thread> renderThreads;

    // 线程同步标志
    std::atomic<bool> isRendering;

public:

    int getCurrentPixelCount() const {
        return metrics.pixelsProcessed.load(std::memory_order_relaxed);
    }

    Zbuffer<float> zbuffer;
    GamesEngineeringBase::Window canvas;

#if defined(USE_SIMD_OPTIMIZATION)
    MatrixSIMD perspective;
#else
    matrix perspective;
#endif

    
    Renderer() {
        canvas.create(1024, 768, "Raster");
        zbuffer.create(1024, 768);
#if defined(USE_SIMD_OPTIMIZATION)
        perspective = MatrixSIMD::makePerspective(fov, aspect, n, f);
#else
        perspective = matrix::makePerspective(fov, aspect, n, f);
#endif
        resetMetrics();
        metrics.activeThreads = 0;
        isRendering = false;
    }
   
    // 添加线程安全的像素记录方法
    void recordPixelProcessed() {
        metrics.pixelsProcessed.fetch_add(1, std::memory_order_relaxed);
    }

    // 修改三角形记录方法也为线程安全
    void recordTriangleProcessed() {
        metrics.triangleCount.fetch_add(1, std::memory_order_relaxed);
    }

    // 添加线程安全的像素绘制函数
    void drawPixelThreadSafe(int x, int y, unsigned char r, unsigned char g, unsigned char b, float depth) {
        std::lock_guard<std::mutex> lock(zbufferMutex);
        if (zbuffer(x, y) > depth) {
            canvas.draw(x, y, r, g, b);
            zbuffer(x, y) = depth;
        }
    }

    // 用于分配渲染任务给线程的结构体
    struct RenderTask {
        int startY;
        int endY;
        int startX;
        int endX;
    };

    // 生成渲染任务的辅助函数
    std::vector<RenderTask> createRenderTasks(int width, int height, int threadCount) {
        std::vector<RenderTask> tasks;
        int rowsPerThread = height / threadCount;

        for (int i = 0; i < threadCount; ++i) {
            RenderTask task;
            task.startY = i * rowsPerThread;
            task.endY = (i == threadCount - 1) ? height : (i + 1) * rowsPerThread;
            task.startX = 0;
            task.endX = width;
            tasks.push_back(task);
        }

        return tasks;
    }

    // 清理线程资源的函数
    void cleanupThreads() {
        isRendering = false;
        for (auto& thread : renderThreads) {
            if (thread.joinable()) {
                thread.join();
            }
        }
        renderThreads.clear();
    }

    // 原有的方法保持不变
    void clear() {
        canvas.clear();
        zbuffer.clear();
        resetMetrics();
    }

    void present() {
        cleanupThreads();  // 确保所有渲染线程完成
        canvas.present();
        updatePerformanceMetrics();
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

    // 性能分析相关函数
    void startFrameMetrics() {
        resetMetrics();
    }

  
    const PerformanceMetrics& getPerformanceMetrics() const {
        return metrics;
    }

    std::string getPerformanceReport() const {
        std::stringstream ss;
        ss << "Performance Report:\n";
        ss << "Frame Time: " << metrics.frameTime << "ms\n";
        ss << "Triangles Processed: " << metrics.triangleCount << "\n";
        ss << "Pixels Processed: " << metrics.pixelsProcessed << "\n";
        ss << "Active Threads: " << metrics.activeThreads << "\n";
        ss << "Average Pixels per Triangle: "
            << (metrics.triangleCount > 0 ? metrics.pixelsProcessed / metrics.triangleCount : 0) << "\n";
        return ss.str();
    }

    

private:
    // 修改重置函数
    void resetMetrics() {
        metrics.reset();
    }

    void updatePerformanceMetrics() {
        static auto lastFrameTime = std::chrono::high_resolution_clock::now();
        auto currentTime = std::chrono::high_resolution_clock::now();
        metrics.frameTime = std::chrono::duration<double, std::milli>(currentTime - lastFrameTime).count();
        lastFrameTime = currentTime;
    }
};