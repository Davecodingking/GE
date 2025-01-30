#pragma once
#include "simd_utils.h"
#include "vec4_simd.h"
#include <cmath>
#include <iostream>
#include <xmmintrin.h> // 用于 SIMD 指令

class alignas(16) MatrixSIMD {
private:
    union {
        __m128 rows[4];           // SIMD友好的行存储
        float elements[16];        // 线性访问数组
    };

public:
    // 默认构造函数：创建单位矩阵
    MatrixSIMD() {
        rows[0] = _mm_set_ps(0.0f, 0.0f, 0.0f, 1.0f);
        rows[1] = _mm_set_ps(0.0f, 0.0f, 1.0f, 0.0f);
        rows[2] = _mm_set_ps(0.0f, 1.0f, 0.0f, 0.0f);
        rows[3] = _mm_set_ps(1.0f, 0.0f, 0.0f, 0.0f);
    }

    // 从16个浮点数构造
    MatrixSIMD(const float* data) {
        for (int i = 0; i < 4; ++i) {
            rows[i] = _mm_load_ps(&data[i * 4]);
        }
    }

    // 矩阵乘法
    MatrixSIMD operator*(const MatrixSIMD& other) const {
        MatrixSIMD result;

        for (int i = 0; i < 4; ++i) {
            __m128 row = rows[i];
            __m128 r0 = _mm_shuffle_ps(row, row, _MM_SHUFFLE(0, 0, 0, 0));
            __m128 r1 = _mm_shuffle_ps(row, row, _MM_SHUFFLE(1, 1, 1, 1));
            __m128 r2 = _mm_shuffle_ps(row, row, _MM_SHUFFLE(2, 2, 2, 2));
            __m128 r3 = _mm_shuffle_ps(row, row, _MM_SHUFFLE(3, 3, 3, 3));

            result.rows[i] = _mm_add_ps(
                _mm_add_ps(
                    _mm_mul_ps(r0, other.rows[0]),
                    _mm_mul_ps(r1, other.rows[1])
                ),
                _mm_add_ps(
                    _mm_mul_ps(r2, other.rows[2]),
                    _mm_mul_ps(r3, other.rows[3])
                )
            );
        }

        return result;
    }

    // 矩阵与向量乘法
    Vec4SIMD operator*(const Vec4SIMD& vec) const {
        __m128 v = vec.getData();
        __m128 result = _mm_setzero_ps();

        for (int i = 0; i < 4; ++i) {
            __m128 row = rows[i];
            __m128 mul = _mm_mul_ps(row, v);
            __m128 sum = _mm_hadd_ps(mul, mul);
            sum = _mm_hadd_ps(sum, sum);

            switch (i) {
            case 0: result = _mm_insert_ps(result, sum, 0 << 4); break;
            case 1: result = _mm_insert_ps(result, sum, 1 << 4); break;
            case 2: result = _mm_insert_ps(result, sum, 2 << 4); break;
            case 3: result = _mm_insert_ps(result, sum, 3 << 4); break;
            }
        }

        return Vec4SIMD(result);
    }


    // 创建平移矩阵
    static MatrixSIMD makeTranslation(float x, float y, float z) {
        MatrixSIMD m;
        m.rows[0] = _mm_set_ps(0.0f, 0.0f, 0.0f, 1.0f);
        m.rows[1] = _mm_set_ps(0.0f, 0.0f, 1.0f, 0.0f);
        m.rows[2] = _mm_set_ps(0.0f, 1.0f, 0.0f, 0.0f);
        m.rows[3] = _mm_set_ps(1.0f, z, y, x);
        return m;
    }

    // 创建X轴旋转矩阵
    static MatrixSIMD makeRotateX(float angleRadians) {
        float c = std::cos(angleRadians);
        float s = std::sin(angleRadians);

        MatrixSIMD m;
        m.rows[0] = _mm_set_ps(0.0f, 0.0f, 0.0f, 1.0f);
        m.rows[1] = _mm_set_ps(0.0f, 0.0f, c, 0.0f);
        m.rows[2] = _mm_set_ps(0.0f, c, 0.0f, 0.0f);
        m.rows[3] = _mm_set_ps(1.0f, s, -s, 0.0f);
        return m;
    }

    // 创建Y轴旋转矩阵
    static MatrixSIMD makeRotateY(float angleRadians) {
        float c = std::cos(angleRadians);
        float s = std::sin(angleRadians);

        MatrixSIMD m;
        m.rows[0] = _mm_set_ps(0.0f, 0.0f, 0.0f, c);
        m.rows[1] = _mm_set_ps(0.0f, 0.0f, 1.0f, 0.0f);
        m.rows[2] = _mm_set_ps(0.0f, c, 0.0f, -s);
        m.rows[3] = _mm_set_ps(1.0f, s, 0.0f, 0.0f);
        return m;
    }

    // 创建Z轴旋转矩阵
    static MatrixSIMD makeRotateZ(float angleRadians) {
        float c = std::cos(angleRadians);
        float s = std::sin(angleRadians);

        MatrixSIMD m;
        m.rows[0] = _mm_set_ps(0.0f, 0.0f, 0.0f, c);
        m.rows[1] = _mm_set_ps(0.0f, 0.0f, c, -s);
        m.rows[2] = _mm_set_ps(0.0f, 1.0f, 0.0f, 0.0f);
        m.rows[3] = _mm_set_ps(1.0f, 0.0f, s, 0.0f);
        return m;
    }

    // 创建缩放矩阵
    static MatrixSIMD makeScale(float s) {
        MatrixSIMD m;
        m.rows[0] = _mm_set_ps(0.0f, 0.0f, 0.0f, s);
        m.rows[1] = _mm_set_ps(0.0f, 0.0f, s, 0.0f);
        m.rows[2] = _mm_set_ps(0.0f, s, 0.0f, 0.0f);
        m.rows[3] = _mm_set_ps(1.0f, 0.0f, 0.0f, 0.0f);
        return m;
    }

    // 创建透视投影矩阵
    static MatrixSIMD makePerspective(float fov, float aspect, float znear, float zfar) {
        float tanHalfFov = std::tan(fov * 0.5f);
        float range = zfar - znear;

        MatrixSIMD m;
        m.rows[0] = _mm_set_ps(0.0f, 0.0f, 0.0f, 1.0f / (aspect * tanHalfFov));
        m.rows[1] = _mm_set_ps(0.0f, 0.0f, 1.0f / tanHalfFov, 0.0f);
        m.rows[2] = _mm_set_ps(0.0f, -(zfar + znear) / range, 0.0f, 0.0f);
        m.rows[3] = _mm_set_ps(-1.0f, -(2.0f * zfar * znear) / range, 0.0f, 0.0f);
        return m;
    }

    // 元素访问操作符
    float& operator()(unsigned int row, unsigned int col) {
        return elements[row * 4 + col];
    }

    // 常量元素访问操作符
    float operator()(unsigned int row, unsigned int col) const {
        return elements[row * 4 + col];
    }

    // 用于调试的显示函数
    void display() const {
        for (int i = 0; i < 4; ++i) {
            alignas(16) float row[4];
            _mm_store_ps(row, rows[i]);
            for (int j = 0; j < 4; ++j) {
                std::cout << row[j] << '\t';
            }
            std::cout << std::endl;
        }
    }
};

#if defined(USE_SIMD_OPTIMIZATION)
using matrix_type = MatrixSIMD;
#else
using matrix_type = matrix;
#endif