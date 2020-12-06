#pragma once

#include <iostream>
#include <string>
#include <string_view>
#include <type_traits>
#include <cassert>
#include <sstream>
#include <vector>
#include <algorithm>
#include <random>
#include <ctime>

template<typename T> std::ostream& operator<<(std::ostream& out, const std::vector<T>& vec) {
	if (vec.empty())
		return out;
	out << vec[0];
	for (std::size_t i = 1; i < vec.size(); ++i) {
		out << " " << vec[i];
	}
	return out;
}
template<typename T> std::istream& operator>>(std::istream&& in, std::vector<T>& vec) {
	std::string line;
	std::getline(in, line);
	std::istringstream ss(line);
	for (T value; ss >> value;) vec.emplace_back(value);
	return in;
}

void write() {}
template<typename... Args> void write(const Args&... args) {
	(std::cout << ... << args);
}
void writef(std::string_view format) {
	write(format);
}
template<typename Arg, typename... Args> void writef(std::string_view format, const Arg& arg, const Args&... rest) {
	std::size_t i = 0;
	for (; i < format.size() && format[i] != '%'; ++i);
	write(format.substr(0, i), arg);
	writef(format.substr(i + 1), rest...);
}
template<typename... Args> void writeln(const Args&... args) {
	write(args...); std::cout << '\n';
}
template<typename... Args> void writefln(std::string_view format, const Args&... args) {
	writef(format, args...); std::cout << '\n';
}

template<typename T=std::string> T read() {
	T value;
	std::cin >> value;
	return value;
}
std::string readln() {
	std::string line;
	std::getline(std::cin, line);
	return line;
}
template<typename T> T parse(std::string_view str) {
	if constexpr (std::is_same_v<T, std::string_view>)
		return str;
	else if constexpr (std::is_same_v<T, std::string>) 
		return std::string{ str };
	else {
		T value;
		std::istringstream(std::string{ str }) >> value;
		return value;
	}
}
template<typename T> std::string toString(const T& value) {
	std::ostringstream ss;
	ss << value;
	return ss.str();
}
template<typename T> T readln() {
	return parse<T>(readln());
}
template<typename T = std::string> std::vector<T> readArray() {
	return readln<std::vector<T>>();
}


//------------
// string utility
//------------
template<typename T=std::string> std::vector<T> split(std::string_view str, std::string_view delim=" ") {
	std::vector<T> result;
	for (std::size_t pos; (pos = str.find(delim)) != std::string::npos;) {
		result.emplace_back(parse<T>(str.substr(0, pos)));
		str = str.substr(pos + delim.size());
	}
	result.emplace_back(parse<T>(str));
	return result;
}
template<typename T> std::string join(const std::vector<T>& vec, std::string sep="") {
	if (vec.size() == 0) 
		return "";
	auto result = toString(vec[0]);
	for (int i = 1; i < vec.size(); ++i) {
		result += sep + toString(vec[i]);
	}
	return result;
}


//------------
// algorithms simplified syntax
//------------
template<typename Container> void sort(Container& container) {
	std::sort(container.begin(), container.end());
}
template<typename Container, typename Pred> void sort(Container& container, Pred pred) {
	std::sort(container.begin(), container.end(), pred);
}


//------------
// Range iterator
//------------
template<typename T> struct SimpleRange {
    struct Iterator {
        Iterator(T value) : value(value) {}
        Iterator& operator++() { ++value; return *this; }
        bool operator==(const Iterator& other) { return value == other.value; }
        bool operator!=(const Iterator& other) { return value != other.value; }
        T operator*() { return value; }
    private:
        T value;
    };

	SimpleRange(T start, T end) : start_(start), end_(end) {}
	SimpleRange(T end) : start_(0), end_(end) {}
    Iterator begin() { return Iterator(start_); }
    Iterator end()   { return Iterator(end_); }
private:
    T start_;
    T end_;
};
template<typename T> struct LinearRange {
    struct Iterator {
        Iterator(T value, T increment) : value(value), increment(increment) {}
        Iterator& operator++() { value += increment; return *this; }
        bool operator==(const Iterator& other) { return value == other.value; }
        bool operator!=(const Iterator& other) { return value != other.value; }
        T operator*() { return value; }
    private:
        T value;
		T increment;
    };

	LinearRange(T start, T end, T increment) : start_(start), end_(end), increment_(increment) {}
    Iterator begin() { return Iterator(start_, increment_); }
    Iterator end()   { return Iterator(end_, increment_); }
private:
    T start_;
    T end_;
	T increment_;
};
template<typename T, typename U> struct GeometricRange {
    struct Iterator {
        Iterator(T value, U multiplyFactor) : value(value), multiplyFactor(multiplyFactor) {}
        Iterator& operator++() { value *= multiplyFactor; return *this; }
        bool operator==(const Iterator& other) { return value == other.value; }
        bool operator!=(const Iterator& other) { return value != other.value; }
        T operator*() { return value; }
    private:
        T value;
		U multiplyFactor;
    };

	GeometricRange(T start, T end, U multiplyFactor) : start_(start), end_(end), multiplyFactor_(multiplyFactor) {}
    Iterator begin() { return Iterator(start_, multiplyFactor_); }
    Iterator end()   { return Iterator(end_, multiplyFactor_); }
private:
    T start_;
    T end_;
	U multiplyFactor_;
};

template<typename Container> auto indicies(const Container& container) -> SimpleRange<decltype(container.size())> {
    return SimpleRange<decltype(container.size())>(0, container.size());
}
template<typename T> auto range(T start, T end, T increment) { return LinearRange(start, end, increment); }
template<typename T> auto range(T start, T end) { return SimpleRange(start, end); }
template<typename T> auto range(T end) { return SimpleRange(end); }
template<typename T, typename U> auto geometricRange(T start, T end, U multiplyFactor) { return GeometricRange(start, end, multiplyFactor); }

//------------
// Random
//------------
static inline std::mt19937_64 _randomGen(std::time(nullptr));
template<typename T = double> T random(T start=0, T end=1) {
	if constexpr (std::is_integral_v<T>) {
		if constexpr (sizeof(T) == 1) {
			std::uniform_int_distribution<int> dist(start, end);
			return static_cast<T>(dist(_randomGen));
		} else {
			std::uniform_int_distribution<T> dist(start, end);
			return dist(_randomGen);
		}
	} else {
		std::uniform_real_distribution<T> dist(start, end);
		return dist(_randomGen);
	}
}
template<typename T = double> std::vector<T> randomVec(uint64_t size) {
	std::vector<T> vec;
	for (auto i : range(size)) {
		if constexpr (std::is_integral_v<T>) vec.emplace_back(random(std::numeric_limits<T>::min(), std::numeric_limits<T>::max()));
		else vec.emplace_back(random());
	}
	return vec;
}



//------------
// 2D array
//------------

template<typename T> struct Array2d {
	std::vector<T> data;
	int rowSize;
	int colSize;

	Array2d(int rowSize, int colSize) : data(rowSize* colSize), rowSize(rowSize), colSize(colSize) {}

	T& operator()(int row, int col) {
		return data[col + row * colSize];
	}
	const T& operator()(int row, int col) const {
		return data[col + row * colSize];
	}
};

