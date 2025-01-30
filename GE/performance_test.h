#pragma once
#include <chrono>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include "renderer.h"
#include <functional>


class PerformanceBenchmark {
private:
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
    };

    std::vector<TestResult> simdResults;
    std::vector<TestResult> baseResults;

public:
    void runBenchmark(const std::string& sceneName,std::function<void()> sceneFunc,
        bool simdEnabled,
        int testDurationSeconds = 10) {
        TestResult result;
        result.sceneName = sceneName;
        result.frameCount = 0;
        result.totalTriangles = 0;
        result.totalPixels = 0;
        result.minFrameTime = (std::numeric_limits<double>::max)();
        result.maxFrameTime = 0;

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

        if (simdEnabled) {
            simdResults.push_back(result);
        }
        else {
            baseResults.push_back(result);
        }

    }

    void generateReport(const std::string& filename) {
        std::ofstream report(filename);
        report << std::fixed << std::setprecision(2);

        report << "Performance Benchmark Report\n";
        report << "===========================\n\n";

        for (size_t i = 0; i < baseResults.size(); i++) {
            const auto& base = baseResults[i];
            const auto& simd = simdResults[i];

            report << "Scene: " << base.sceneName << "\n";
            report << "----------------------------------------\n";
            report << "                   Base        SIMD     Improvement\n";
            report << "Average FPS:       " << std::setw(8) << base.averageFPS
                << "    " << std::setw(8) << simd.averageFPS
                << "    " << std::setw(8) << (simd.averageFPS / base.averageFPS * 100 - 100) << "%\n";

            report << "Frame Time (ms):   " << std::setw(8) << base.averageFrameTime
                << "    " << std::setw(8) << simd.averageFrameTime
                << "    " << std::setw(8) << (base.averageFrameTime / simd.averageFrameTime * 100 - 100) << "%\n";

            report << "Min Frame Time:    " << std::setw(8) << base.minFrameTime
                << "    " << std::setw(8) << simd.minFrameTime << "\n";
            report << "Max Frame Time:    " << std::setw(8) << base.maxFrameTime
                << "    " << std::setw(8) << simd.maxFrameTime << "\n";

            report << "\n";
        }
    }
};