#pragma once
#include <chrono>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include "renderer.h"
#include <functional>

// 定义优化类型枚举
enum OptimizationType {
    BASE,           // 基础版本
    SIMD,           // SIMD优化版本
    THREADED,       // 多线程优化版本
    SIMD_THREADED   // SIMD+多线程组合优化
};

class PerformanceBenchmark {
public:
    struct TestResult {
        std::string sceneName;
        double averageFrameTime;
        double minFrameTime;
        double maxFrameTime;
        double totalTime;
        int frameCount;
        int totalTriangles;
        int totalPixels;
        double averageFPS;
        OptimizationType type;  // 添加优化类型字段

        // 新增字段
        int trianglesPerFrame;    // 每帧平均三角形数
        int pixelsPerFrame;       // 每帧平均像素数
        double peakPerformance;   // 峰值性能（最高FPS）
  
        double pixelsPerMs;       // 每毫秒处理的像素数
       
    };

    // 存储所有类型的测试结果
    std::vector<TestResult> results;

public:
    void runBenchmark(const std::string& sceneName,
        std::function<void()> sceneFunc,
        OptimizationType type,
        int testDurationSeconds = 10) {
        TestResult result;
        result.sceneName = sceneName;
        result.frameCount = 0;
        result.totalTriangles = 0;
        result.totalPixels = 0;
        result.minFrameTime = (std::numeric_limits<double>::max)();
        result.maxFrameTime = 0;
        result.type = type;

        auto startTime = std::chrono::high_resolution_clock::now();
        auto endTime = startTime + std::chrono::seconds(testDurationSeconds);

        while (std::chrono::high_resolution_clock::now() < endTime) {
            auto frameStart = std::chrono::high_resolution_clock::now();

            

            sceneFunc();

            auto frameEnd = std::chrono::high_resolution_clock::now();
            double frameTime = std::chrono::duration<double, std::milli>(frameEnd - frameStart).count();

            result.minFrameTime = (std::min)(result.minFrameTime, frameTime);
            result.maxFrameTime = (std::max)(result.maxFrameTime, frameTime);
            result.frameCount++;
        }

        auto testEnd = std::chrono::high_resolution_clock::now();
        result.totalTime = std::chrono::duration<double, std::milli>(testEnd - startTime).count();
        result.averageFrameTime = result.totalTime / result.frameCount;
        result.averageFPS = 1000.0 / result.averageFrameTime;

        results.push_back(result);
    }

    void generateReport(const std::string& filename) {
        std::ofstream report(filename);
        report << std::fixed << std::setprecision(2);

        report << "Performance Benchmark Report\n";
        report << "===========================\n\n";

        // 按场景分组显示结果
        std::string currentScene = "";
        std::vector<TestResult> sceneResults;

        for (const auto& result : results) {
            if (currentScene != result.sceneName) {
                if (!sceneResults.empty()) {
                    printSceneComparison(report, sceneResults);
                    sceneResults.clear();
                }
                currentScene = result.sceneName;
            }
            sceneResults.push_back(result);
        }

        // 打印最后一个场景的结果
        if (!sceneResults.empty()) {
            printSceneComparison(report, sceneResults);
        }
    }

private:
    void printSceneComparison(std::ofstream& report, const std::vector<TestResult>& sceneResults) {
    report << "Scene: " << sceneResults[0].sceneName << "\n";
    report << "----------------------------------------\n";
    report << std::left;
    
    // 表头
    report << std::setw(12) << "Version"
           << std::setw(10) << "FPS"
           << std::setw(12) << "Frame Time"
           << std::setw(15) << "Triangles/f"    // 新增：每帧三角形数
           << std::setw(12) << "Pixels/f"       // 新增：每帧像素数
           << std::setw(15) << "Pixels/ms"      // 新增：像素处理速率
           << std::setw(12) << "Improvement"
           << "\n";

    // 找到基础版本用于计算改进百分比
    auto baseResult = std::find_if(sceneResults.begin(), sceneResults.end(),
        [](const TestResult& r) { return r.type == BASE; });

    for (const auto& result : sceneResults) {
        report << std::setw(12);

        // 版本名称
        switch (result.type) {
        case BASE: report << "Base"; break;
        case SIMD: report << "SIMD"; break;
        case THREADED: report << "Threaded"; break;
        case SIMD_THREADED: report << "SIMD+Thread"; break;
        }

        // 基础性能指标
        report << std::fixed << std::setprecision(2)
               << std::setw(10) << result.averageFPS
               << std::setw(12) << result.averageFrameTime;

        // 新增的详细性能指标
        report << std::setw(15) << result.trianglesPerFrame
               << std::setw(12) << result.pixelsPerFrame
               << std::setw(15) << (result.pixelsPerFrame / result.averageFrameTime);

        // 性能改进百分比
        if (baseResult != sceneResults.end() && result.type != BASE) {
            double improvement = (baseResult->averageFrameTime / result.averageFrameTime - 1.0) * 100;
            report << std::setw(12) << "+" << improvement << "%";
        }
        report << "\n";
    }

    // 添加场景总结
    report << "\nScene Summary:\n";
    report << "--------------\n";
    report << "Total Test Duration: " << std::fixed << std::setprecision(2) 
           << sceneResults[0].totalTime / 1000.0 << " seconds\n";
    report << "Average Scene Complexity: " << sceneResults[0].totalTriangles << " triangles\n";
    report << "Total Pixels Processed: " << sceneResults[0].totalPixels << "\n\n";

    // 添加性能分析
    report << "Performance Analysis:\n";
    report << "--------------------\n";
    if (sceneResults.size() > 1) {
        auto& base = *baseResult;
        auto& best = *std::max_element(sceneResults.begin(), sceneResults.end(),
            [](const TestResult& a, const TestResult& b) {
                return a.averageFPS < b.averageFPS;
            });
        
        report << "Best Performance: " << best.averageFPS << " FPS ("
               << (best.type == BASE ? "Base" :
                   best.type == SIMD ? "SIMD" :
                   best.type == THREADED ? "Threaded" : "SIMD+Thread")
               << " version)\n";
        report << "Overall Improvement: " 
               << ((best.averageFPS / base.averageFPS - 1.0) * 100) 
               << "% over base version\n";
        /*report << "Pixel Processing Rate: " 
               << (best.pixelsPerFrame / best.averageFrameTime) 
               << " pixels/ms (peak)\n\n";*/
    }
}
};