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
        report << "Version      FPS    Frame Time   Improvement\n";

        // 找到基础版本的结果用于计算改进百分比
        auto baseResult = std::find_if(sceneResults.begin(), sceneResults.end(),
            [](const TestResult& r) { return r.type == BASE; });

        for (const auto& result : sceneResults) {
            report << std::left << std::setw(12);

            // 输出版本名称
            switch (result.type) {
            case BASE: report << "Base"; break;
            case SIMD: report << "SIMD"; break;
            case THREADED: report << "Threaded"; break;
            case SIMD_THREADED: report << "SIMD+Thread"; break;
            }

            report << std::fixed << std::setprecision(2)
                << std::setw(8) << result.averageFPS
                << std::setw(12) << result.averageFrameTime;

            // 计算改进百分比
            if (baseResult != sceneResults.end() && result.type != BASE) {
                double improvement = (baseResult->averageFrameTime / result.averageFrameTime - 1.0) * 100;
                report << std::setw(8) << "+" << improvement << "%";
            }
            report << "\n";
        }
        report << "\n";
    }
};