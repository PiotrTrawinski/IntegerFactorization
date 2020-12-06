#pragma once
#include <immintrin.h>
#include <type_traits>
#include <new>
#include <cstdint>
#include <limits>
#include "compilerMacros.h"
#include <iostream>

namespace avx {
    template<typename T> constexpr bool IsCompatible = 
        ((sizeof(T) == 1 || sizeof(T) == 2 || sizeof(T) == 4 || sizeof(T) == 8) && std::is_integral_v<T>) 
        || std::numeric_limits<T>::is_iec559;

#ifdef __AVX2__
    #define AVX2_IS_AVAILABLE
    constexpr bool IsAvailable = true;        
#else
    constexpr bool IsAvailable = false;
#endif
}

template<typename ValueType, typename Enable = void> struct AvxType;

#ifdef AVX2_IS_AVAILABLE

// 8, 16, 32 or 64 bit integer type
template<typename ValueType> struct AvxType<ValueType, typename std::enable_if_t<(sizeof(ValueType) == 1 || sizeof(ValueType) == 2 || sizeof(ValueType) == 4 || sizeof(ValueType) == 8) && std::is_integral_v<ValueType>>> {
    __m256i value;
    AvxType(const __m256i& value) : value(value) {}
    operator __m256i() const { return value; }
};

// IEEE 754 binary32 type
template<> struct AvxType<float, typename std::enable_if_t<sizeof(float) == 4 && std::numeric_limits<float>::is_iec559>> {
    __m256 value;
    AvxType(const __m256& value) : value(value) {}
    operator __m256() const { return value; }
};

// IEEE 754 binary64 type
template<> struct AvxType<double, typename std::enable_if_t<sizeof(double) == 8 && std::numeric_limits<double>::is_iec559>> {
    __m256d value;
    AvxType(const __m256d& value) : value(value) {}
    operator __m256d() const { return value; }
};

#else
template<> struct AvxType<float> {};
template<> struct AvxType<double> {};
#endif



namespace avx {
    template<typename T> constexpr int packedCount() {
        return 256 / 8 / sizeof(T);
    }

#ifdef AVX2_IS_AVAILABLE

    // load/store from/to memory
    template<typename T> AvxType<T> loadAligned(const T* ptr) { return _mm256_load_si256((__m256i*)ptr); }
    AvxType<float> loadAligned(const float* ptr)              { return _mm256_load_ps(ptr); }
    AvxType<double> loadAligned(const double* ptr)            { return _mm256_load_pd(ptr); }

    template<typename T> AvxType<T> loadUnaligned(const T* ptr) { return _mm256_loadu_si256((__m256i*)ptr); }
    AvxType<float> loadUnaligned(const float* ptr)              { return _mm256_loadu_ps(ptr); }
    AvxType<double> loadUnaligned(const double* ptr)            { return _mm256_loadu_pd(ptr); }

    template<typename T> AvxType<T> load(const T* ptr) { return avx::loadUnaligned(ptr); }

    template<typename T> void storeAligned(T* dst, const AvxType<T>& src) { return _mm256_store_si256((__m256i*)dst, src); }
    void storeAligned(float* dst, AvxType<float> src)              { return _mm256_store_ps(dst, src); }
    void storeAligned(double* dst, AvxType<double> src)            { return _mm256_store_pd(dst, src); }

    template<typename T> void storeUnaligned(T* dst, const AvxType<T>& src) { return _mm256_storeu_si256((__m256i*)dst, src); }
    void storeUnaligned(float* dst, AvxType<float> src)              { return _mm256_storeu_ps(dst, src); }
    void storeUnaligned(double* dst, AvxType<double> src)            { return _mm256_storeu_pd(dst, src); }


    // addition
    template<typename T> AvxType<T> add(const AvxType<T>& a, const AvxType<T>& b) {
        if constexpr (sizeof(T) == 1) return _mm256_add_epi8(a, b);
        if constexpr (sizeof(T) == 2) return _mm256_add_epi16(a, b);
        if constexpr (sizeof(T) == 4) return _mm256_add_epi32(a, b);
        if constexpr (sizeof(T) == 8) return _mm256_add_epi64(a, b);
        return __m256i{}; // impossible to hit, but some compilerers can't see that so needed to avoid warnings 
    }
    AvxType<float> add(const AvxType<float>& a, const AvxType<float>& b)    { return _mm256_add_ps(a, b); }
    AvxType<double> add(const AvxType<double>& a, const AvxType<double>& b) { return _mm256_add_pd(a, b); }

    // subtraction
    template<typename T> AvxType<T> sub(const AvxType<T>& a, const AvxType<T>& b) {
        if constexpr (sizeof(T) == 1) return _mm256_sub_epi8(a, b);
        if constexpr (sizeof(T) == 2) return _mm256_sub_epi16(a, b);
        if constexpr (sizeof(T) == 4) return _mm256_sub_epi32(a, b);
        if constexpr (sizeof(T) == 8) return _mm256_sub_epi64(a, b);
        return __m256i{}; // impossible to hit, but some compilerers can't see that so needed to avoid warnings 
    }
    AvxType<float> sub(const AvxType<float>& a, const AvxType<float>& b)    { return _mm256_sub_ps(a, b); }
    AvxType<double> sub(const AvxType<double>& a, const AvxType<double>& b) { return _mm256_sub_pd(a, b); }

    // multiplication
    template<typename T> AvxType<T> mul(const AvxType<T>& a, const AvxType<T>& b) {
        if constexpr (sizeof(T) == 1) {
            auto even = _mm256_mullo_epi16(a, b);
            auto odd = _mm256_mullo_epi16(_mm256_srli_epi16(a, 8), _mm256_srli_epi16(b, 8));
            return _mm256_or_si256(_mm256_slli_epi16(odd, 8), _mm256_and_si256(even, _mm256_set1_epi16(0xFF)));
        }
        if constexpr (sizeof(T) == 2) return _mm256_mullo_epi16(a, b);
        if constexpr (sizeof(T) == 4) return _mm256_mullo_epi32(a, b);
        if constexpr (sizeof(T) == 8) {
            // taken from https://stackoverflow.com/a/37320416
            __m256i bswap = _mm256_shuffle_epi32(b, 0xB1);
            __m256i prodlh = _mm256_mullo_epi32(a, bswap);

            __m256i prodlh2 = _mm256_srli_epi64(prodlh, 32);
            __m256i prodlh3 = _mm256_add_epi32(prodlh2, prodlh);
            __m256i prodlh4 = _mm256_and_si256(prodlh3, _mm256_set1_epi64x(0x00000000FFFFFFFF));

            __m256i prodll = _mm256_mul_epu32(a, b);
            __m256i prod = _mm256_add_epi64(prodll, prodlh4);
            return  prod;
        }
        return __m256i{}; // impossible to hit, but some compilerers can't see that so needed to avoid warnings 
    }
    AvxType<float> mul(const AvxType<float>& a, const AvxType<float>& b)    { return _mm256_mul_ps(a, b); }
    AvxType<double> mul(const AvxType<double>& a, const AvxType<double>& b) { return _mm256_mul_pd(a, b); }

    // division
    AvxType<float> div(const AvxType<float>& a, const AvxType<float>& b) { return _mm256_div_ps(a, b); }
    AvxType<double> div(const AvxType<double>& a, const AvxType<double>& b) { return _mm256_div_pd(a, b); }


    // initialization
    template<typename T> AvxType<T> zero() { return _mm256_setzero_si256(); }
    template<> AvxType<float> zero()       { return _mm256_setzero_ps(); }
    template<> AvxType<double> zero()      { return _mm256_setzero_pd(); }

    template<typename T> AvxType<T> setElements(T value) {
        if constexpr (sizeof(T) == 1) return _mm256_set1_epi8(value);
        if constexpr (sizeof(T) == 2) return _mm256_set1_epi16(value);
        if constexpr (sizeof(T) == 4) return _mm256_set1_epi32(value);
        if constexpr (sizeof(T) == 8) return _mm256_set1_epi64x(value);
        return __m256i{}; // impossible to hit, but some compilerers can't see that so needed to avoid warnings 
    }
    template<> AvxType<float> setElements(float value)   { return _mm256_set1_ps(value); }
    template<> AvxType<double> setElements(double value) { return _mm256_set1_pd(value); }

    template<typename T> AvxType<T> setElements(T v1, T v2, T v3, T v4, T v5, T v6, T v7, T v8, T v9, T v10, T v11, T v12, T v13, T v14, T v15, T v16, T v17, T v18, T v19, T v20, T v21, T v22, T v23, T v24, T v25, T v26, T v27, T v28, T v29, T v30, T v31, T v32) {
        return _mm256_setr_epi8(v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16, v17, v18, v19, v20, v21, v22, v23, v24, v25, v26, v27, v28, v29, v30, v31, v32);
    }
    template<typename T> AvxType<T> setElements(T v1, T v2, T v3, T v4, T v5, T v6, T v7, T v8, T v9, T v10, T v11, T v12, T v13, T v14, T v15, T v16) {
        return _mm256_setr_epi16(v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16);
    }
    template<typename T> AvxType<T> setElements(T v1, T v2, T v3, T v4, T v5, T v6, T v7, T v8) {
        return _mm256_setr_epi32(v1, v2, v3, v4, v5, v6, v7, v8);
    }
    template<typename T> AvxType<T> setElements(T v1, T v2, T v3, T v4) {
        return _mm256_setr_epi64x(v1, v2, v3, v4);
    }
    template<> AvxType<float> setElements(float v1, float v2, float v3, float v4, float v5, float v6, float v7, float v8) {
        return _mm256_setr_ps(v1, v2, v3, v4, v5, v6, v7, v8); 
    }
    template<> AvxType<double> setElements(double v1, double v2, double v3, double v4) {
        return _mm256_setr_pd(v1, v2, v3, v4);
    }

    // get single elements
    template<int Index, typename T> T getElement(const AvxType<T>& a) {
    #if defined(COMPILER_GCC) || defined(COMPILER_CLANG)
        return a.value[Index];
    #else
        if constexpr (sizeof(T) == 8) return _mm256_extract_epi64(a, Index);
        if constexpr (sizeof(T) == 4) return _mm256_extract_epi32(a, Index);
        if constexpr (sizeof(T) == 2) return _mm256_extract_epi16(a, Index);
        if constexpr (sizeof(T) == 1) return _mm256_extract_epi8(a, Index);
        return 0;
    #endif
    }
    template<int Index> float getElement(const AvxType<float>& a) {
    #if defined(COMPILER_GCC) || defined(COMPILER_CLANG)
        return a.value[Index];
    #else
        return 0; 
    #endif
    }
    template<int Index> double getElement(const AvxType<double>& a) {
    #if defined(COMPILER_GCC) || defined(COMPILER_CLANG)
        return a.value[Index];
    #else
        switch (Index) {
        case 0: return _mm256_cvtsd_f64(a);
        case 1: return _mm256_cvtsd_f64(_mm256_unpackhi_pd(a, a));
        case 2: return _mm_cvtsd_f64(_mm256_extractf128_pd(a, 1));
        case 3: { auto b = _mm256_extractf128_pd(a, 1); return _mm_cvtsd_f64(_mm_unpackhi_pd(b, b)); }
        default: return 0;
        }
    #endif
    }


    // other
    AvxType<float> floor(const AvxType<float>& a) {
        return _mm256_floor_ps(a);;
    }
    AvxType<double> floor(const AvxType<double>& a) {
        return _mm256_floor_pd(a);;
    }
    AvxType<double> equal(const AvxType<double>& a, const AvxType<double>& b) {
        return _mm256_cmp_pd(a, b, _CMP_EQ_OQ);
    }
    AvxType<double> greaterOrEqual(const AvxType<double>& a, const AvxType<double>& b) {
        return _mm256_cmp_pd(a, b, _CMP_GE_OQ);
    }
    AvxType<double> lessOrEqual(const AvxType<double>& a, const AvxType<double>& b) {
        return _mm256_cmp_pd(a, b, _CMP_LE_OQ);
    }
    template<typename T> uint32_t mask8bit(const AvxType<T>& a) {
        if constexpr (std::is_same_v<T, double>) return _mm256_movemask_epi8(_mm256_castpd_si256(a));
        else if constexpr (std::is_same_v<T, float>) return _mm256_movemask_epi8(_mm256_castps_si256(a));
        else return _mm256_movemask_epi8(a);
    }
#else
    template<typename T> AvxType<T> loadAligned(const T* ptr) { return AvxType<T>{}; }
    template<typename T> AvxType<T> loadUnaligned(const T* ptr) { return AvxType<T>{}; }
    template<typename T> void storeAligned(T* dst, const AvxType<T>& src) {}
    template<typename T> void storeUnaligned(T* dst, const AvxType<T>& src) {}
    template<typename T> AvxType<T> add(const AvxType<T>& a, const AvxType<T>& b) { return AvxType<T>{}; }
    template<typename T> AvxType<T> sub(const AvxType<T>& a, const AvxType<T>& b) { return AvxType<T>{}; }
    template<typename T> AvxType<T> mul(const AvxType<T>& a, const AvxType<T>& b) { return AvxType<T>{}; }
    template<typename T> AvxType<T> div(const AvxType<T>& a, const AvxType<T>& b) { return AvxType<T>{}; }
    template<typename T> AvxType<T> zero() { return AvxType<T>{}; }
    template<typename T> AvxType<T> setElements(T value) { return AvxType<T>{}; }
#endif
}

template<typename T> AvxType<T> operator+(const AvxType<T>& a, const AvxType<T>& b) {
    return avx::add(a, b);
}
template<typename T> AvxType<T> operator-(const AvxType<T>& a, const AvxType<T>& b) {
    return avx::sub(a, b);
}
template<typename T> AvxType<T> operator*(const AvxType<T>& a, const AvxType<T>& b) {
    return avx::mul(a, b);
}
template<typename T> AvxType<T> operator/(const AvxType<T>& a, const AvxType<T>& b) {
    return avx::div(a, b);
}
template<typename T> AvxType<T> operator+=(const AvxType<T>& a, const AvxType<T>& b) {
    return a = a + b;
}
template<typename T> AvxType<T> operator-=(const AvxType<T>& a, const AvxType<T>& b) {
    return a = a - b;
}
template<typename T> AvxType<T> operator*=(const AvxType<T>& a, const AvxType<T>& b) {
    return a = a * b;
}
template<typename T> AvxType<T> operator/=(const AvxType<T>& a, const AvxType<T>& b) {
    return a = a / b;
}
template<typename T> AvxType<T> operator==(const AvxType<T>& a, const AvxType<T>& b) {
    return avx::equal(a, b);
}
template<typename T> AvxType<T> operator<=(const AvxType<T>& a, const AvxType<T>& b) {
    return avx::lessOrEqual(a, b);
}
template<typename T> AvxType<T> operator>=(const AvxType<T>& a, const AvxType<T>& b) {
    return avx::greaterOrEqual(a, b);
}
