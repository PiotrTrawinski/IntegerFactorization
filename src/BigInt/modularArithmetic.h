#pragma once

#include <cassert>
#include "../Utility/alwaysInline.h"
#include "../Utility/debugAssert.h"
#include "BigIntMaxCap.h"
#include "BigIntFixedSize.h"
#include "BigIntGmp.h"
#include "Unsigned64.h"

namespace modArithm {
    /*
        requires value type T and modular type U to have defined free functions:
        void assign(T& dest, const T& src)
        void assign(T& dest, const int& src)
        void modAdd(T& result, const T& a, const T& b, const U& mod)
        void modSub(T& result, const T& a, const T& b, const U& mod)
        void modMul(T& result, const T& a, const T& b, const U& mod)
        void modSqr(T& result, const T& a, const U& mod)
        void modNeg(T& result, const T& a, const U& mod)
        void modInv(T& result, const T& a, const U& mod)
        void modDbl(T& result, const T& a, const U& mod)
        void shr(T& result, const T& a, int c)
    */

template <typename T> struct Expression {
    ALWAYS_INLINE const auto& mod()   const      { return static_cast<const T&>(*this).mod(); }
    ALWAYS_INLINE const auto& value() const      { return static_cast<const T&>(*this).value(); }
    ALWAYS_INLINE auto& value()                  { return static_cast<T&>(*this).value(); }
    template<typename V> ALWAYS_INLINE bool allOperandsDifferent(Expression<V>& ref) const         { return static_cast<const T&>(*this).allOperandsDifferent(ref); }
    template<typename V> ALWAYS_INLINE bool allOperandsButFirstDifferent(Expression<V>& ref) const { return static_cast<const T&>(*this).allOperandsButFirstDifferent(ref); }
    template<typename E> ALWAYS_INLINE void eval(Expression<E>& result) const { static_cast<const T&>(*this).eval(result); }
};

template<typename T, typename U> struct Value : public Expression<Value<T, U>> {
    constexpr static bool IsValue = true;

    ALWAYS_INLINE Value(T& value, const U& mod) : value_(value), mod_(mod) {}

    ALWAYS_INLINE const U& mod()   const { return mod_; }
    ALWAYS_INLINE const T& value() const { return value_; }
    ALWAYS_INLINE T& value() { return value_; }
    template<typename V> ALWAYS_INLINE bool allOperandsDifferent(Expression<V>& ref) const         { return &value_ != &ref.value(); }
    template<typename V> ALWAYS_INLINE bool allOperandsButFirstDifferent(Expression<V>& ref) const { return true; }

    template<typename E> ALWAYS_INLINE void eval(Expression<E>& result) const {
        assign(result.value(), value());
    }
    template<typename E> ALWAYS_INLINE Value& operator=(const Expression<E>& exp) {
        exp.eval(*this);
        return *this;
    }
    ALWAYS_INLINE Value& operator=(const int& val) {
        assign(value(), val);
        return *this;
    }
    ALWAYS_INLINE Value& operator=(const Value& exp) {
        exp.eval(*this);
        return *this;
    }
    ALWAYS_INLINE Value& operator=(const T& val) {
        assign(value(), val);
        return *this;
    }
    ALWAYS_INLINE Value& operator>>=(const int& val) {
        shr(value(), value(), val);
        return *this;
    }
private:
    T& value_;
    const U& mod_;
};
template<typename T, typename U> struct ConstValue : public Expression<ConstValue<T, U>> {
    constexpr static bool IsValue = true;

    ALWAYS_INLINE ConstValue(const T& value, const U& mod) : value_(value), mod_(mod) {}
    
    ALWAYS_INLINE const U& mod()   const { return mod_; }
    ALWAYS_INLINE const T& value() const { return value_; }

    template<typename V> ALWAYS_INLINE bool allOperandsDifferent(Expression<V>& ref) const { return &value_ != &ref.value(); }
    template<typename V> ALWAYS_INLINE bool allOperandsButFirstDifferent(Expression<V>& ref) const { return true; }
    
    template<typename E> ALWAYS_INLINE void eval(Expression<E>& result) const {
        assign(result.value(), value());
    }
private:
    const T& value_;
    const U& mod_;
};


enum class Operation {
    Add, Sub, Mul, Sqr, Neg, Inv, Dbl
};
template <Operation Op, typename T1, typename T2> struct BinaryOperation : public Expression<BinaryOperation<Op, T1, T2>> {
    constexpr static bool IsValue = false;

    ALWAYS_INLINE BinaryOperation(const Expression<T1>& a, const Expression<T2>& b) : a(a), b(b) {
        debugAssert(&a.mod() == &b.mod());
    }

    ALWAYS_INLINE const auto& mod() const {
        return a.mod();
    }

    template<typename V> ALWAYS_INLINE bool allOperandsDifferent(Expression<V>& ref) const { 
        return a.allOperandsDifferent(ref) && b.allOperandsDifferent(ref);
    }
    template<typename V> ALWAYS_INLINE bool allOperandsButFirstDifferent(Expression<V>& ref) const { 
        return a.allOperandsButFirstDifferent(ref) && b.allOperandsDifferent(ref);
    }

    template<typename T> ALWAYS_INLINE void eval(Expression<T>& result) const {
        static_assert(T2::IsValue, "right side of expression has to be a simple value (not expression)");
        debugAssert(allOperandsButFirstDifferent(result), "result value can only be used as first operand in first operation (otherwise silent temporary variables would need to be made)");
        if constexpr (!T1::IsValue) {
            a.eval(result);
            if constexpr (Op == Operation::Add) modAdd(result.value(), result.value(), b.value(), result.mod());
            if constexpr (Op == Operation::Sub) modSub(result.value(), result.value(), b.value(), result.mod());
            if constexpr (Op == Operation::Mul) modMul(result.value(), result.value(), b.value(), result.mod());
        } else {
            if constexpr (Op == Operation::Add) modAdd(result.value(), a.value(), b.value(), result.mod());
            if constexpr (Op == Operation::Sub) modSub(result.value(), a.value(), b.value(), result.mod());
            if constexpr (Op == Operation::Mul) modMul(result.value(), a.value(), b.value(), result.mod());
        }
    }

private:
    const Expression<T1>& a;
    const Expression<T2>& b;
};
template <Operation Op, typename T1> struct UnaryOperation : public Expression<UnaryOperation<Op, T1>> {
    constexpr static bool IsValue = false;

    ALWAYS_INLINE UnaryOperation(const Expression<T1>& a) : a(a) {}

    ALWAYS_INLINE const auto& mod() const {
        return a.mod();
    }

    template<typename V> ALWAYS_INLINE bool allOperandsDifferent(Expression<V>& ref) const {
        return a.allOperandsDifferent(ref);
    }
    template<typename V> ALWAYS_INLINE bool allOperandsButFirstDifferent(Expression<V>& ref) const {
        return a.allOperandsButFirstDifferent(ref);
    }

    template<typename T> ALWAYS_INLINE void eval(Expression<T>& result) const {
        debugAssert(allOperandsButFirstDifferent(result), "result value can only be used as first operand in first operation (otherwise silent temporary variables would need to be made)");
        if constexpr (!T1::IsValue) {
            a.eval(result);
            if constexpr (Op == Operation::Sqr) modSqr(result.value(), result.value(), result.mod());
            if constexpr (Op == Operation::Neg) modNeg(result.value(), result.value(), result.mod());
            if constexpr (Op == Operation::Inv) modInv(result.value(), result.value(), result.mod());
            if constexpr (Op == Operation::Dbl) modDbl(result.value(), result.value(), result.mod());
        } else {
            if constexpr (Op == Operation::Sqr) modSqr(result.value(), a.value(), result.mod());
            if constexpr (Op == Operation::Neg) modNeg(result.value(), a.value(), result.mod());
            if constexpr (Op == Operation::Inv) modInv(result.value(), a.value(), result.mod());
            if constexpr (Op == Operation::Dbl) modDbl(result.value(), a.value(), result.mod());
        }
    }

private:
    const Expression<T1>& a;
};
template <typename T1, typename T2> BinaryOperation<Operation::Add, T1, T2> ALWAYS_INLINE operator+(const Expression<T1>& a, const Expression<T2>& b) {
    return BinaryOperation<Operation::Add, T1, T2>(a, b);
}
template <typename T1, typename T2> BinaryOperation<Operation::Sub, T1, T2> ALWAYS_INLINE operator-(const Expression<T1>& a, const Expression<T2>& b) {
    return BinaryOperation<Operation::Sub, T1, T2>(a, b);
}
template <typename T1, typename T2> BinaryOperation<Operation::Mul, T1, T2> ALWAYS_INLINE operator*(const Expression<T1>& a, const Expression<T2>& b) {
    return BinaryOperation<Operation::Mul, T1, T2>(a, b);
}
template <typename T1> UnaryOperation<Operation::Sqr, T1> ALWAYS_INLINE sqr(const Expression<T1>& a) {
    return UnaryOperation<Operation::Sqr, T1>(a);
}
template <typename T1> UnaryOperation<Operation::Neg, T1> ALWAYS_INLINE operator-(const Expression<T1>& a) {
    return UnaryOperation<Operation::Neg, T1>(a);
}
template <typename T1> UnaryOperation<Operation::Inv, T1> ALWAYS_INLINE inv(const Expression<T1>& a) {
    return UnaryOperation<Operation::Inv, T1>(a);
}
template <typename T1> UnaryOperation<Operation::Dbl, T1> ALWAYS_INLINE dbl(const Expression<T1>& a) {
    return UnaryOperation<Operation::Dbl, T1>(a);
}

template<typename T, typename... Args> inline auto createContext(const T& mod, Args&... args) {
    return std::make_tuple(modArithm::Value(args, mod)...);
}

}
