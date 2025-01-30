#pragma once
#include "simd_utils.h"
#include "vec4.h"
#include <iostream>

// 确保类按16字节对齐以优化SIMD操作
class alignas(16) Vec4SIMD {
private:
    __m128 data;  // 使用SSE 4浮点数向量存储数据

public:
    // 默认构造函数：初始化为(0,0,0,1)
    Vec4SIMD() : data(_mm_set_ps(1.0f, 0.0f, 0.0f, 0.0f)) {}

    // 从分量构造
    Vec4SIMD(float x, float y, float z, float w = 1.0f)
        : data(_mm_set_ps(w, z, y, x)) {
    }

    // 直接从__m128构造
    explicit Vec4SIMD(__m128 m) : data(m) {}

    // 从vec4构造
    explicit Vec4SIMD(const vec4& v)
        : data(_mm_set_ps(v[3], v[2], v[1], v[0])) {
    }

    // 获取内部SIMD数据
    __m128 getData() const { return data; }

    // SIMD优化的点乘操作
    static float dot(const Vec4SIMD& v1, const Vec4SIMD& v2) {
        __m128 mul = _mm_mul_ps(v1.data, v2.data);
        __m128 hadd1 = _mm_hadd_ps(mul, mul);
        __m128 hadd2 = _mm_hadd_ps(hadd1, hadd1);
        return _mm_cvtss_f32(hadd2);
    }

    // SIMD优化的向量标准化
    void normalise() {
        // 计算平方和
        __m128 square = _mm_mul_ps(data, data);
        __m128 sum = _mm_hadd_ps(square, square);
        sum = _mm_hadd_ps(sum, sum);

        // 使用rsqrt指令计算倒数平方根（快速但精度较低）
        __m128 invLength = _mm_rsqrt_ps(sum);

        // Newton-Raphson迭代提高精度
        // y = y * (3 - x * y * y) / 2
        __m128 half = _mm_set1_ps(0.5f);
        __m128 three = _mm_set1_ps(3.0f);
        __m128 y2 = _mm_mul_ps(invLength, invLength);
        __m128 xy2 = _mm_mul_ps(sum, y2);
        __m128 three_minus_xy2 = _mm_sub_ps(three, xy2);
        invLength = _mm_mul_ps(_mm_mul_ps(invLength, three_minus_xy2), half);

        // 归一化向量
        data = _mm_mul_ps(data, invLength);

        // 确保w分量为1
        data = _mm_blend_ps(data, _mm_set_ps(1.0f, 0.0f, 0.0f, 0.0f), 0x8);
    }

    // SIMD优化的叉乘
    Vec4SIMD cross(const Vec4SIMD& other) const {
        __m128 a = data;
        __m128 b = other.data;

        __m128 a_yzx = _mm_shuffle_ps(a, a, _MM_SHUFFLE(3, 0, 2, 1));
        __m128 b_yzx = _mm_shuffle_ps(b, b, _MM_SHUFFLE(3, 0, 2, 1));
        __m128 a_zxy = _mm_shuffle_ps(a, a, _MM_SHUFFLE(3, 1, 0, 2));
        __m128 b_zxy = _mm_shuffle_ps(b, b, _MM_SHUFFLE(3, 1, 0, 2));

        __m128 result = _mm_sub_ps(
            _mm_mul_ps(a_yzx, b_zxy),
            _mm_mul_ps(a_zxy, b_yzx)
        );

        // 设置w分量为0（叉积结果应该是向量）
        result = _mm_blend_ps(result, _mm_setzero_ps(), 0x8);

        return Vec4SIMD(result);
    }

    // 除以w分量
    void divideW() {
        __m128 w = _mm_shuffle_ps(data, data, _MM_SHUFFLE(3, 3, 3, 3));
        data = _mm_div_ps(data, w);
        // 重置w分量为1
        data = _mm_blend_ps(data, _mm_set_ps(1.0f, 0.0f, 0.0f, 0.0f), 0x8);
    }

    // 向量加法
    Vec4SIMD operator+(const Vec4SIMD& other) const {
        return Vec4SIMD(_mm_add_ps(data, other.data));
    }

    // 向量减法
    Vec4SIMD operator-(const Vec4SIMD& other) const {
        return Vec4SIMD(_mm_sub_ps(data, other.data));
    }

    // 标量乘法
    Vec4SIMD operator*(float scalar) const {
        return Vec4SIMD(_mm_mul_ps(data, _mm_set1_ps(scalar)));
    }

    // 分量访问
    float& operator[](const unsigned int index) {
        alignas(16) float temp[4];
        _mm_store_ps(temp, data);
        return temp[index];
    }

    float operator[](const unsigned int index) const {
        alignas(16) float temp[4];
        _mm_store_ps(temp, data);
        return temp[index];
    }

    // 转换回vec4
    operator vec4() const {
        alignas(16) float temp[4];
        _mm_store_ps(temp, data);
        return vec4(temp[0], temp[1], temp[2], temp[3]);
    }

    // 调试用显示函数
    void display() const {
        alignas(16) float temp[4];
        _mm_store_ps(temp, data);
        std::cout << temp[0] << '\t' << temp[1] << '\t'
            << temp[2] << '\t' << temp[3] << std::endl;
    }
};

// 如果需要全局启用SIMD优化，可以使用类型别名
#if defined(USE_SIMD_OPTIMIZATION)
using vec4_type = Vec4SIMD;
#else
using vec4_type = vec4;
#endif