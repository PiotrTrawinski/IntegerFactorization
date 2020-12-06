#pragma once
#include "../Utility/compilerMacros.h"
#include "../Utility/debugAssert.h"
#include "../Utility/alwaysInline.h"
#include "64bitIntrinsics.h"
#include <cstdint>
#include <string>
#include <type_traits>

/*
    Wrapper classes for simple bit set/get operations.
*/
struct BigIntBit {
    uint64_t* ptr;
    uint32_t index;

    BigIntBit(uint64_t* ptr, uint32_t index) : ptr(ptr), index(index) {}
    void setBit(bool x) {
        *ptr ^= (-x ^ *ptr) & (1ull << index);
    }
    BigIntBit& operator=(bool x) {
        setBit(x);
        return *this;
    }
    BigIntBit& operator=(const BigIntBit& other) {
        setBit(other);
        return *this;
    }
    operator bool() const {
        return (*ptr >> index) & 1;
    }
};
struct BigIntConstBit {
    const uint64_t* ptr;
    uint32_t index;
    BigIntConstBit(const uint64_t* ptr, uint32_t index) : ptr(ptr), index(index) {}
    operator bool() const {
        return (*ptr >> index) & 1;
    }
};
BigIntBit bit(uint64_t& value, uint32_t index) {
    return BigIntBit(&value, index);
}
BigIntConstBit bit(const uint64_t& value, uint32_t index) {
    return BigIntConstBit(&value, index);
}

/*
    Base for string parsing. All big integer types should specialize this
*/
template<typename T, typename = void> struct BigIntParseImpl {
    static T parse(const std::string& str) {
        throw std::logic_error("method is not implemented for given type");
    }
};
template<typename T> T bigIntParse(const std::string& str) {
    return BigIntParseImpl<T>::parse(str);
}


/* 
    Representation of value, which allows for faster modulo operation (by means of multiplication + shift).
    It's known as 'Barrett reduction', see: 
        Barrett, P. (1986). 
        "Implementing the Rivest Shamir and Adleman Public Key Encryption Algorithm on a Standard Digital Signal Processor". 
        Advances in Cryptology — CRYPTO' 86. Lecture Notes in Computer Science. 263. pp. 311–323
*/
template<typename T> struct BarretReductionMod {
    T mod;           // the original value
    T R;             // R = 2^k / mod
    mutable T q_tmp; // temporary for speed/convienience (temporary value is required in order to perform modulo operation)
    uint32_t k;      // k = log2(mod) * 2
};


template<typename T> struct MontgomeryReductionMod {
    MontgomeryReductionMod() {}
    MontgomeryReductionMod(const T& mod) : mod(mod) {}
    T mod;           // the original value
    T k;             // k = (r*(r^-1 % mod) - 1) / mod, where r = 2^b
    uint32_t b;      // b = ceilToLimb(log2(mod))
    mutable T s;     // temporary used in calculations
    mutable T t;     // temporary used in calculations
};

template<typename Test, template<typename...> class Ref> struct is_specialization : std::false_type {};
template<template<typename...> class Ref, typename... Args> struct is_specialization<Ref<Args...>, Ref> : std::true_type {};

template<typename T> using BigIntValueType = 
    typename std::conditional<is_specialization<T, BarretReductionMod>::value, decltype(T::mod), 
    typename std::conditional<is_specialization<T, MontgomeryReductionMod>::value, decltype(T::mod),
    T
>::type>::type;
