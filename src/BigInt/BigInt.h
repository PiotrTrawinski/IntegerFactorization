#pragma once
#include <variant>
#include "BigIntFixedSize.h"
#include "BigIntGmp.h"
#include "Unsigned64.h"

struct BigInt {
    BigInt() {}
    BigInt(uint64_t n) {
        _montgomeryNumberData = MontgomeryReductionMod<BigIntFixedSize<1>>(n);
    }
    BigInt(const BigIntGmp& n) {
        _montgomeryNumberData = MontgomeryReductionMod<BigIntGmp>(n);
    }
    template<int Size> BigInt(const BigIntFixedSize<Size>& n) {
        _montgomeryNumberData = MontgomeryReductionMod<BigIntFixedSize<Size>>(n);
    }
    BigInt(BigIntGmp&& n) {
        _montgomeryNumberData = MontgomeryReductionMod<BigIntGmp>(std::move(n));
    }
    template<int Size> BigInt(BigIntFixedSize<Size>&& n) {
        _montgomeryNumberData = MontgomeryReductionMod<BigIntFixedSize<Size>>(std::move(n));
    }
    BigInt(const MontgomeryReductionMod<BigIntGmp>& n) {
        _montgomeryNumberData = n;
    }
    template<int Size> BigInt(const MontgomeryReductionMod<BigIntFixedSize<Size>>& n) {
        _montgomeryNumberData = n;
    }
    BigInt(MontgomeryReductionMod<BigIntGmp>&& n) {
        _montgomeryNumberData = std::move(n);
    }
    template<int Size> BigInt(MontgomeryReductionMod<BigIntFixedSize<Size>>&& n) {
        _montgomeryNumberData = std::move(n);
    }
    BigInt(const std::string& str) {
        initializeModValue(str);
    }
    BigInt& operator=(const BigInt& n) {
        _montgomeryNumberData = n._montgomeryNumberData;
        isMontgomeryFormInitialized = n.isMontgomeryFormInitialized;
        return *this;
    }



    template<typename Lambda> auto visit(Lambda lambda) const {
        initializeIfNeeded();
        return std::visit(lambda, _montgomeryNumberData);
    }
    template<typename Lambda> auto visit(Lambda lambda) {
        initializeIfNeeded();
        return std::visit(lambda, _montgomeryNumberData);
    }
    template<typename Lambda> auto visitNoInit(Lambda lambda) const {
        return std::visit(lambda, _montgomeryNumberData);
    }
    template<typename Lambda> auto visitNoInit(Lambda lambda) {
        return std::visit(lambda, _montgomeryNumberData);
    }

    std::string toString() const {
        return visit([](auto&& n) { return ::toString(n.mod); });
    }

    template<typename T, typename U> static constexpr bool IsType() {
        return std::is_same_v<std::remove_cv_t<std::remove_reference_t<T>>, MontgomeryReductionMod<U>>;
    }

    template<typename T> bool IsType() const {
        return std::holds_alternative<MontgomeryReductionMod<T>>(_montgomeryNumberData);
    }

    template<typename T> T& get() {
        return std::get<MontgomeryReductionMod<T>>(_montgomeryNumberData).mod;
    }
    template<typename T> const T& get() const {
        return std::get<MontgomeryReductionMod<T>>(_montgomeryNumberData).mod;
    }

    int sizeInBits() const {
        return visitNoInit([](auto&& a) { return a.mod.sizeInBits(); });
    }

    int typeIndex() const {
        return _montgomeryNumberData.index();
    }

    bool isOne() const {
        return visitNoInit([](auto&& a) { return a.mod == decltype(a.mod){ 1 }; });
    }
    
    void fitToSize() {
        auto expectedIndex = expectedTypeIndex(sizeInBits());
        auto index = typeIndex();
        if (index == expectedIndex) {
            return;
        }
        switch (expectedIndex) {
        case 0: fitToSize<BigIntGmp>(); break;
        case 1: fitToSize<BigIntFixedSize<1>>(); break;
        case 2: fitToSize<BigIntFixedSize<2>>(); break;
        case 3: fitToSize<BigIntFixedSize<3>>(); break;
        case 4: fitToSize<BigIntFixedSize<4>>(); break;
        case 5: fitToSize<BigIntFixedSize<5>>(); break;
        case 6: fitToSize<BigIntFixedSize<6>>(); break;
        case 7: fitToSize<BigIntFixedSize<7>>(); break;
        case 8: fitToSize<BigIntFixedSize<8>>(); break;
        }
    }

    BigInt& operator/=(const BigInt& n) {
        visitNoInit([&n](auto&& a) { n.visitNoInit([&a](auto&& b) {
            div(a.mod, a.mod, decltype(a.mod)(b.mod));
            //a /= b;
        });});
        isMontgomeryFormInitialized = false;
        fitToSize();
        return *this;
    }

private:
    int expectedTypeIndex(int sizeInBits) {
        if (sizeInBits < 64)  return 1;
        if (sizeInBits < 128) return 2;
        if (sizeInBits < 192) return 3;
        if (sizeInBits < 256) return 4;
        if (sizeInBits < 320) return 5;
        if (sizeInBits < 384) return 6;
        if (sizeInBits < 448) return 7;
        if (sizeInBits < 512) return 8;
        else return 0;
    }

    void initializeIfNeeded() const {
        if (!isMontgomeryFormInitialized) {
            visitNoInit([&](auto&& a) {
                _montgomeryNumberData = getMontgomeryReductionMod(a.mod);
            });
        }
    }
    template<int S> void initializeModValueFixed(const std::string& str) {
        _montgomeryNumberData = MontgomeryReductionMod<BigIntFixedSize<S>>{};
        std::get<MontgomeryReductionMod<BigIntFixedSize<S>>>(_montgomeryNumberData).mod = BigIntFixedSize<S>{ str };
    }
    void initializeModValue(const std::string& str) {
        if      (str.size() <=  18) initializeModValueFixed<1>(str);
        else if (str.size() <=  38) initializeModValueFixed<2>(str);
        else if (str.size() <=  57) initializeModValueFixed<3>(str);
        else if (str.size() <=  76) initializeModValueFixed<4>(str);
        else if (str.size() <=  96) initializeModValueFixed<5>(str);
        else if (str.size() <= 115) initializeModValueFixed<6>(str);
        else if (str.size() <= 134) initializeModValueFixed<7>(str);
        else if (str.size() <= 153) initializeModValueFixed<8>(str);
        else {
            _montgomeryNumberData = MontgomeryReductionMod<BigIntGmp>{};
            std::get<MontgomeryReductionMod<BigIntGmp>>(_montgomeryNumberData).mod = BigIntGmp{ str };
        }
    }
    template<typename T> void fitToSize() {
        visitNoInit([&](auto&& n) {
            _montgomeryNumberData = MontgomeryReductionMod<T>(T(n.mod));
        });
        isMontgomeryFormInitialized = false;
    }

    mutable std::variant<
        MontgomeryReductionMod<BigIntGmp>,
        MontgomeryReductionMod<BigIntFixedSize<1>>,
        MontgomeryReductionMod<BigIntFixedSize<2>>,
        MontgomeryReductionMod<BigIntFixedSize<3>>,
        MontgomeryReductionMod<BigIntFixedSize<4>>,
        MontgomeryReductionMod<BigIntFixedSize<5>>,
        MontgomeryReductionMod<BigIntFixedSize<6>>,
        MontgomeryReductionMod<BigIntFixedSize<7>>,
        MontgomeryReductionMod<BigIntFixedSize<8>>
    > _montgomeryNumberData;
    bool isMontgomeryFormInitialized = false;
};

bool operator==(const BigInt& a, const BigInt& b) {
    if (a.typeIndex() != b.typeIndex())
        return false;
    switch (a.typeIndex()) {
    case 0: return a.get<BigIntGmp>() == b.get<BigIntGmp>();
    case 1: return a.get<BigIntFixedSize<1>>() == b.get<BigIntFixedSize<1>>();
    case 2: return a.get<BigIntFixedSize<2>>() == b.get<BigIntFixedSize<2>>();
    case 3: return a.get<BigIntFixedSize<3>>() == b.get<BigIntFixedSize<3>>();
    case 4: return a.get<BigIntFixedSize<4>>() == b.get<BigIntFixedSize<4>>();
    case 5: return a.get<BigIntFixedSize<5>>() == b.get<BigIntFixedSize<5>>();
    case 6: return a.get<BigIntFixedSize<6>>() == b.get<BigIntFixedSize<6>>();
    case 7: return a.get<BigIntFixedSize<7>>() == b.get<BigIntFixedSize<7>>();
    case 8: return a.get<BigIntFixedSize<8>>() == b.get<BigIntFixedSize<8>>();
    }
    return false;
}
bool operator!=(const BigInt& a, const BigInt& b) {
    return !(a == b);
}

std::ostream& operator<<(std::ostream& out, const BigInt& a) {
    return out << a.toString();
}

BigInt operator ""_big(const char* str) {
    return BigInt(std::string(str));
}
