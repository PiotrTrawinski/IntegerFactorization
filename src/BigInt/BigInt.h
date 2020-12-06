#pragma once
#include <variant>
#include "BigIntFixedSize.h"
#include "BigIntGmp.h"
#include "Unsigned64.h"

struct BigInt {
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

private:
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

std::ostream& operator<<(std::ostream& out, const BigInt& a) {
    return out << a.toString();
}

BigInt operator ""_big(const char* str) {
    return BigInt(std::string(str));
}
