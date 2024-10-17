#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <cmath>
#include <concepts>
#include <string>

#include "data-block.hpp"
#include "matrix.hpp"
#include "operation-constraints.hpp"

// for annotating functions
struct Preserve { };
using Return = Preserve;
struct Modify { };
using NoReturn = Modify;

template <typename T, size_t N> 
struct RawVector;

template <typename T, size_t N>
std::string ToString(RawVector<T, N> const& crVec) {
	return ToString(crVec.mData);
}

template <typename T, size_t N>
struct RawVector {
private:
	DataBlock<T, N> mData;
public:
	template <std::convertible_to<T> U>
	explicit RawVector(DataBlock<U, N> const& crDB) : mData(crDB) { }

	RawVector(std::convertible_to<T> auto const&... args) requires (sizeof...(args) == N) : mData(args...)  { }

	RawVector(T const& crArg) {
		for (int i = 0; i < N; ++i) {
			mData[i] = crArg;
		}
	}

	template <typename U> requires (CanSum<T, U>)
	inline RawVector<sum_t<T, U>, N> operator+(RawVector<U, N> const& crRHS) {
		return RawVector<sum_t<T, U>, N>(mData + crRHS.mData);
	}

	template <typename U> requires (CanSubtract<T, U>)
	inline RawVector<difference_t<T, U>, N> operator-(RawVector<U, N> const& crRHS) {
		return RawVector<difference_t<T, U>, N>(mData - crRHS.mData);
	}

	template <typename U> requires (CanMultiply<T, U>)
	inline RawVector<product_t<T, U>, N> operator*(RawVector<U, N> const& crRHS) {
		return RawVector<product_t<T, U>, N>(mData * crRHS.mData);
	}

	template <typename U> requires (CanDivide<T, U>)
	inline RawVector<quotient_t<T, U>, N> operator/(RawVector<U, N> const& crRHS) {
		return RawVector<quotient_t<T, U>, N>(mData / crRHS.mData);
	}

	template <typename U> requires (CanModulus<T, U>)
	inline RawVector<modulus_t<T, U>, N> operator%(RawVector<U, N> const& crRHS) {
		return RawVector<modulus_t<T, U>, N>(mData % crRHS.mData);
	}

	template <typename U> requires (CanBitwiseXor<T, U>)
	inline RawVector<bitwise_xor_t<T, U>, N> operator^(RawVector<U, N> const& crRHS) {
		return RawVector<bitwise_xor_t<T, U>, N>(mData ^ crRHS.mData);
	}

	template <typename U> requires (CanBitwiseAnd<T, U>)
	inline RawVector<bitwise_and_t<T, U>, N> operator&(RawVector<U, N> const& crRHS) {
		return RawVector<bitwise_and_t<T, U>, N>(mData & crRHS.mData);
	}

	template <typename U> requires (CanBitwiseOr<T, U>)
	inline RawVector<bitwise_or_t<T, U>, N> operator|(RawVector<U, N> const& crRHS) {
		return RawVector<bitwise_or_t<T, U>, N>(mData | crRHS.mData);
	}

	template <typename U = T> requires (CanBitwiseNot<U>)
	inline RawVector<bitwise_not_t<U>, N> operator~() {
		return RawVector<bitwise_not_t<U>, N>(~mData);
	}

	template <std::convertible_to<T> U>
	inline RawVector<T, N>& operator=(RawVector<U, N> const& crRHS) {
		mData = crRHS.mData;
		return *this;
	}

	template <typename U> 
	requires (std::convertible_to<sum_t<T, U>, T>)
	inline RawVector<T, N>& operator+=(RawVector<U, N> const& crRHS) {
		mData += crRHS.mData;

		return *this;
	}

	template <typename U>
	requires (std::convertible_to<difference_t<T, U>, T>)
	inline RawVector<T, N>& operator-=(RawVector<U, N> const& crRHS) {
		mData -= crRHS.mData;

		return *this;
	}

	template <typename U>
	requires (std::convertible_to<product_t<T, U>, T>)
	inline RawVector<T, N>& operator*=(RawVector<U, N> const& crRHS) {
		mData *= crRHS.mData;

		return *this;
	}

	template <typename U>
	requires (std::convertible_to<quotient_t<T, U>, T>)
	inline RawVector<T, N>& operator/=(RawVector<U, N> const& crRHS) {
		mData /= crRHS.mData;

		return *this;
	}

	template <typename U>
	requires (std::convertible_to<modulus_t<T, U>, T>)
	inline RawVector<T, N>& operator%=(RawVector<U, N> const& crRHS) {
		mData %= crRHS.mData;

		return *this;
	}

	template <typename U>
	requires (std::convertible_to<bitwise_xor_t<T, U>, T>)
	inline RawVector<T, N>& operator^=(RawVector<U, N> const& crRHS) {
		mData ^= crRHS.mData;

		return *this;
	}

	template <typename U>
	requires (std::convertible_to<bitwise_and_t<T, U>, T>)
	inline RawVector<T, N>& operator&=(RawVector<U, N> const& crRHS) {
		mData &= crRHS.mData;

		return *this;
	}

	template <typename U>
	requires (std::convertible_to<bitwise_or_t<T, U>, T>)
	inline RawVector<T, N>& operator|=(RawVector<U, N> const& crRHS) {
		mData |= crRHS.mData;

		return *this;
	}

	template <typename U> requires (CanLeftShift<T, U>)
	inline RawVector<left_shift_t<T, U>, N> operator<<(RawVector<U, N> const& crRHS) {
		return RawVector<left_shift_t<T, U>, N>(mData << crRHS.mData);
	}

	template <typename U> requires (CanRightShift<T, U>)
	inline RawVector<right_shift_t<T, U>, N> operator>>(RawVector<U, N> const& crRHS) {
		return RawVector<right_shift_t<T, U>, N>(mData >> crRHS.mData);
	}

	template <typename U> requires (CanLogicalEquality<T, U>)
	inline bool operator==(RawVector<U, N> const& crRHS) {
		return mData == crRHS.mData;
	}

	template <typename U> requires (CanLogicalInequality<T, U>)
	inline bool operator!=(RawVector<U, N> const& crRHS) {
		return mData != crRHS.mData;
	}

	template <CanPrefixIncrement U = T>
	inline RawVector<U, N>& operator++() {
		++mData;

		return *this;
	}

	template <CanPostfixIncrement U = T>
	inline RawVector<U, N> operator++(int) {
		RawVector<U, N> copy = *this;

		mData++;

		return copy;
	}

	template <CanPrefixDecrement U = T>
	inline RawVector<U, N>& operator--() {
		--mData;

		return *this;
	}

	template <CanPostfixDecrement U = T>
	inline RawVector<U, N> operator--(int) {
		RawVector<U, N> copy = *this;
		mData--;

		return copy;
	}

	inline T& operator()(size_t i) {
		return mData[i];
	}

	inline T const& operator()(size_t i) const {
		return mData[i];
	}

	inline T& operator[](size_t i) {
		return mData[i];
	}

	inline T const& operator[](size_t i) const {
		return mData[i];
	}

	friend std::string ToString<T, N>(RawVector<T, N> const& crVec);
};

template <typename T, size_t N>
struct AnnotatedVector : public RawVector<T, N> {

private:
	int mfFlags;
public:

};

template <typename T, size_t N>
using Vector = DataBlock<T, N>;

template <size_t N>
using bvec = Vector<bool, N>;

using bvec2 = Vector<bool, 2>;
using bvec3 = Vector<bool, 3>;
using bvec4 = Vector<bool, 4>;

template <size_t N>
using ivec = Vector<int, N>;

using ivec2 = Vector<int, 2>;
using ivec3 = Vector<int, 3>;
using ivec4 = Vector<int, 4>;

template <size_t N>
using uvec = Vector<unsigned int, N>;

using uvec2 = Vector<unsigned int, 2>;
using uvec3 = Vector<unsigned int, 3>;
using uvec4 = Vector<unsigned int, 4>;

template <size_t N>
using vec = Vector<float, N>;

using vec2 = Vector<float, 2>;
using vec3 = Vector<float, 3>;
using vec4 = Vector<float, 4>;

template <size_t N>
using dvec = Vector<double, N>;

using dvec2 = Vector<double, 2>;
using dvec3 = Vector<double, 3>;
using dvec4 = Vector<double, 4>;

template <typename T, size_t N>
using VectorView = DataBlockView<T, N>;

template <size_t N>
using bvec_view = VectorView<bool, N>;

using bvec2_view = VectorView<bool, 2>;
using bvec3_view = VectorView<bool, 3>;
using bvec4_view = VectorView<bool, 4>;

template <size_t N>
using ivec_view = VectorView<int, N>;

using ivec2_view = VectorView<int, 2>;
using ivec3_view = VectorView<int, 3>;
using ivec4_view = VectorView<int, 4>;

template <size_t N>
using uvec_view = VectorView<unsigned int, N>;

using uvec2_view = VectorView<unsigned int, 2>;
using uvec3_view = VectorView<unsigned int, 3>;
using uvec4_view = VectorView<unsigned int, 4>;

template <size_t N>
using vec_view = VectorView<float, N>;

using vec2_view = VectorView<float, 2>;
using vec3_view = VectorView<float, 3>;
using vec4_view = VectorView<float, 4>;

template <size_t N>
using dvec_view = VectorView<double, N>;

using dvec2_view = VectorView<double, 2>;
using dvec3_view = VectorView<double, 3>;
using dvec4_view = VectorView<double, 4>;

template <typename T, size_t N>
using ConstVectorView = DataBlockView<const T, N>;

template <size_t N>
using const_bvec_view = ConstVectorView<bool, N>;

using const_bvec2_view = ConstVectorView<bool, 2>;
using const_bvec3_view = ConstVectorView<bool, 3>;
using const_bvec4_view = ConstVectorView<bool, 4>;

template <size_t N>
using const_ivec_view = ConstVectorView<int, N>;

using const_ivec2_view = ConstVectorView<int, 2>;
using const_ivec3_view = ConstVectorView<int, 3>;
using const_ivec4_view = ConstVectorView<int, 4>;

template <size_t N>
using const_uvec_view = ConstVectorView<unsigned int, N>;

using const_uvec2_view = ConstVectorView<unsigned int, 2>;
using const_uvec3_view = ConstVectorView<unsigned int, 3>;
using const_uvec4_view = ConstVectorView<unsigned int, 4>;

template <size_t N>
using const_vec_view = ConstVectorView<float, N>;

using const_vec2_view = ConstVectorView<float, 2>;
using const_vec3_view = ConstVectorView<float, 3>;
using const_vec4_view = ConstVectorView<float, 4>;

template <size_t N>
using const_dvec_view = ConstVectorView<double, N>;

using const_dvec2_view = ConstVectorView<double, 2>;
using const_dvec3_view = ConstVectorView<double, 3>;
using const_dvec4_view = ConstVectorView<double, 4>;

template <typename T>
concept Numeric = std::is_arithmetic_v<T>;

template <typename T, size_t N = 0>
inline constexpr bool IsVectorTypeV = false;

template <typename T, size_t N0>
inline constexpr bool IsVectorTypeV<Vector<T, N0>, 0> = true;

template <typename T, size_t N0>
inline constexpr bool IsVectorTypeV<VectorView<T, N0>, 0> = true;

template <typename T, size_t N>
inline constexpr bool IsVectorTypeV<Vector<T, N>, N> = true;

template <typename T, size_t N>
inline constexpr bool IsVectorTypeV<VectorView<T, N>, N> = true;

template <typename T, size_t N = 0>
concept VectorType = IsVectorTypeV<T, N>;

//template <typename T, size_t R, size_t C>
//inline constexpr bool IsMatrixTypeV = false;
//
//template <typename T, size_t R, size_t C>
//inline constexpr bool IsMatrixTypeV<Matrix<T, R, C>, R, C> = true;
//
//template <typename T, size_t R, size_t C>
//inline constexpr bool IsMatrixTypeV<MatrixView<T, R, C>, R, C> = true;
//
//template <typename T, size_t R, size_t C>
//concept MatrixType = IsMatrixTypeV<T, R, C>;

template <typename T>
using ValueType = typename std::remove_cvref_t<T>::value_type;


// might change to a field requirement later
template <Numeric T, size_t N>
T SquaredMagnitude(Vector<T, N> const& crVec) {
	T result{};

	for (int i = 0; i < N; ++i) {
		result += crVec[i] * crVec[i];
	}

	return result;
}

template <typename T>
concept ReturnDescriptor = (std::is_same_v<T, Return> || std::is_same_v<T, NoReturn>);

template <Numeric T, size_t N>
T Magnitude(Vector<T, N> const& crVec) {
	return sqrt(SquaredMagnitude(crVec));
}

template <ReturnDescriptor RD = Preserve>
auto Normalize(VectorType auto& rVec) -> std::conditional_t<std::is_same_v<RD, NoReturn>, void, std::remove_cvref_t<decltype(rVec)>> {
	using T = ValueType<decltype(rVec)>;

	T magnitude_inv = static_cast<T>(1) / Magnitude(rVec);

	if constexpr (std::is_same_v<RD, NoReturn>) {
		rVec *= magnitude_inv;
	}

	else {
		return rVec * magnitude_inv;
	}
}

auto Normalize(VectorType auto const& crVec) -> std::remove_cvref_t<decltype(crVec)> {
	using T = ValueType<decltype(crVec)>;

	T magnitude_inv = static_cast<T>(1) / Magnitude(crVec);

	return crVec * magnitude_inv;
}


template <Numeric T, size_t N>
T Dot(Vector<T, N> const& crLHS, Vector<T, N> const& crRHS) {
	T result{};

	for (int i = 0; i < N; ++i) {
		result += crLHS[i] * crRHS[i];
	}

	return result;
}

template <Numeric T, size_t N>
T Dot(VectorView<T, N> const& crLHS, VectorView<T, N> const& crRHS) {
	std::remove_cvref_t<T> result{};

	for (int i = 0; i < N; ++i) {
		result += crLHS[i] * crRHS[i];
	}

	return result;
}

template <Numeric T>
Vector<T, 3> Cross(Vector<T, 3> const& crLHS, Vector<T, 3> const& crRHS) {
	Vector<T, 3> result;

	result[0] = crLHS[1] * crRHS[2] - crRHS[1] * crLHS[2];
	result[1] = crRHS[0] * crLHS[2] - crLHS[0] * crRHS[2];
	result[2] = crLHS[0] * crRHS[1] - crRHS[0] * crLHS[1];

	return result;
}

template <Numeric T>
Vector<std::remove_cvref_t<T>, 3> Cross(VectorView<T, 3> const& crLHS, VectorView<T, 3> const& crRHS) {
	Vector<std::remove_cvref_t<T>, 3> result;

	result[0] = crLHS[1] * crRHS[2] - crRHS[1] * crLHS[2];
	result[1] = crRHS[0] * crLHS[2] - crLHS[0] * crRHS[2];
	result[2] = crLHS[0] * crRHS[1] - crRHS[0] * crLHS[1];

	return result;
}

template <Numeric T, size_t N>
bool AreOrthogonal(Vector<T, N> const& crVec) {
	return true;
}

template <Numeric T, Numeric... Ts, size_t N>
bool AreOrthogonal(Vector<T, N> const& crVec, Vector<Ts, N> const&... crVecs) {
	return ((Dot(crVec, crVecs) == 0) && ...) && AreOrthogonal(crVecs...);
}

template <Numeric T, size_t N>
bool IsUnit(Vector<T, N> const& crVec) {
	return SquaredMagnitude(crVec) == 1;
}

template <Numeric... Ts, size_t N>
bool AreOrthonormal(Vector<Ts, N> const &... crVecs) {
	return AreOrthogonal(crVecs...) && (IsUnit(crVecs) && ...);
}
#endif