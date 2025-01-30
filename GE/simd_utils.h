#pragma once
#include <immintrin.h>
#include <intrin.h>

class SIMDUtils {
public:
    // 检查CPU是否支持必要的SIMD指令集
    static bool checkSIMDSupport() {
        int cpuInfo[4];
        __cpuid(cpuInfo, 1);

        // 检查SSE、AVX支持
        bool hasSSE = cpuInfo[3] & (1 << 25);
        bool hasAVX = cpuInfo[2] & (1 << 28);

        // 检查AVX2支持
        __cpuid(cpuInfo, 7);
        bool hasAVX2 = cpuInfo[1] & (1 << 5);

        return hasSSE && hasAVX && hasAVX2;
    }

    // 内存对齐辅助函数
    static void* alignedAlloc(size_t size) {
        return _aligned_malloc(size, 32); // AVX需要32字节对齐
    }

    static void alignedFree(void* ptr) {
        _aligned_free(ptr);
    }
};

// 定义对齐宏
#define SIMD_ALIGN alignas(32)