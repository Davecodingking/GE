#pragma once
#include "vec4_test.h"

class Vec4OperationsTest : public Vec4Test {
public:
    // 测试点乘运算的正确性
    bool testDotProduct() override {
        std::cout << "\nTesting Dot Product:\n";

        // 使用原始版本计算点乘
        float originalResult = vec4::dot(testData.original1, testData.original2);

        // 使用SIMD版本计算点乘
        float simdResult = vec4_simd::dot(testData.simd1, testData.simd2);

        // 比较结果
        bool success = isEqual(originalResult, simdResult);

        // 输出测试结果
        std::cout << "Original result: " << originalResult << "\n";
        std::cout << "SIMD result: " << simdResult << "\n";
        std::cout << "Dot Product Test: " << (success ? "PASSED" : "FAILED") << "\n";

        return success;
    }

    // 测试点乘运算的性能
    void benchmarkDotProduct() override {
        std::cout << "\nBenchmarking Dot Product:\n";

        const int iterations = 1000000;  // 百万次运算以获得可靠的性能数据

        // 测试原始版本性能
        double originalTime = measureTime([&]() {
            volatile float result = vec4::dot(testData.original1, testData.original2);
            }, iterations);

        // 测试SIMD版本性能
        double simdTime = measureTime([&]() {
            volatile float result = vec4_simd::dot(testData.simd1, testData.simd2);
            }, iterations);

        // 输出性能比较结果
        std::cout << "Original version took: " << originalTime << " ms\n";
        std::cout << "SIMD version took: " << simdTime << " ms\n";
        std::cout << "Speedup: " << originalTime / simdTime << "x\n";
    }

    // 可以继续添加其他运算的测试...
};