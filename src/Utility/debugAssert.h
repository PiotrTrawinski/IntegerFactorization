#pragma once
#include <iostream>
#include <string>
#include <stdexcept>

namespace details {
    void printErrorValues() {}

    template<typename T, typename... Ts> 
    void printErrorValues(const T& value, const Ts&... values) {
        std::cerr << value;
        printErrorValues(values...);
    }

    template<typename... Ts>
    void assertInfo(const char* expression, const char* filePath, int lineNumber, const Ts&... values) {
        std::cerr << "Assertion failed\n";
        std::cerr << "file       : " << filePath << '\n';
        std::cerr << "line       : " << lineNumber << '\n';
        std::cerr << "expression : " << expression << '\n';
        if (sizeof...(Ts) > 0) {
            std::cerr << "message    : ";
            printErrorValues(values...);
            std::cerr << '\n';
        }
        std::cerr << std::flush;
    }
}

#ifdef _DEBUG

#define debugAssert(expression, ...)\
do {\
    if (!(expression)) {\
        ::details::assertInfo(#expression, __FILE__, __LINE__, __VA_ARGS__);\
        throw std::logic_error(#expression);\
    }\
} while (false)

#else

#define debugAssert(expression, ...)\
do {\
    (void)sizeof(expression);\
} while (false)

#endif
