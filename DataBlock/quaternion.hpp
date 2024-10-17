//#ifndef QUATERNION_HPP
//#define QUATERNION_HPP
//
//#include "data-block.hpp"
//#include "vector.hpp"
//
//namespace Quat {
//
//	template <Numeric T>
//	
//
//	template <Numeric T>
//	Vector<T, 4> Mul(Vector<T, 4> const& crLHS, Vector<T, 4> const& crRHS) {
//		ConstVectorView<T, 3> imLHS = QuatIm(crLHS);
//		ConstVectorView<T, 3> imRHS = QuatIm(crRHS);
//
//		T const& wLHS = Re(crLHS);
//		T const& wRHS = Re(crRHS);
//
//		return { wLHS * wRHS - Dot(imLHS, imRHS), Cross(imLHS, imRHS) + wLHS * imRHS + wRHS * imLHS };
//	}
//}
//#endif

#ifndef QUATERNION_HPP
#define QUATERNION_HPP

#include <cstddef>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>

#include "data-block.hpp"
#include "index-range.hpp"
#include "matrix.hpp"
#include "vector.hpp"

template <Numeric T>
struct Quaternion;

template <Numeric T>
Quaternion<T> Mul(Quaternion<T> const& crLHS, Quaternion<T> const& crRHS);

template <Numeric T>
struct Quaternion {
private:
	typedef T Quaternion::* const _q[4];
	static _q const q;

public: 
	T w, x, y, z;

	T& operator()(size_t i) {
		return this->*q[i];
	}

	T const& operator()(size_t i) const {
		return this->*q[i];
	}

	T& operator[](size_t i) {
		return this->*q[i]; 
	}

	T const& operator[](size_t i) const {
		return this->*q[i];
	}

	Quaternion<T>& operator=(Quaternion<T> const& crRHS) {
		memcpy(this, &crRHS, sizeof(Quaternion<T>));
		return *this;
	}

	Quaternion(T _w, Vector<T, 3> const& crVec) : w{ _w }, x{ crVec[0] }, y{ crVec[1] },
		z{ crVec[2] } { }

	Quaternion(T _w, VectorView<T, 3> const& crVecView) : w{ _w }, x{ crVecView[0] }, 
		y{ crVecView[1] }, z{ crVecView[2] } { }

	Quaternion<T> operator*(Quaternion<T> const& crRHS) {
		return Mul(*this, crRHS);
	}

	Quaternion<T>& operator*=(Quaternion<T> const& crRHS) {
		*this = Mul(*this, crRHS); 
		return *this;
	}
};

template <Numeric T>
Quaternion<T>::_q const Quaternion<T>::q = { &Quaternion<T>::w,
&Quaternion<T>::x, &Quaternion<T>::y, &Quaternion<T>::z };

template <Numeric T>
std::string ToString(Quaternion<T> const& crQuat, int width = 3) {
	std::ostringstream stream;

	stream << "Re: " << Width(Re(crQuat), width) << '\n';
	stream << "Im: " << ToString(Im(crQuat), width);

	return stream.str();
}

template <Numeric T>
T& Re(Quaternion<T>& rQuat) {
	return rQuat.w;
}

template <Numeric T>
T const& Re(Quaternion<T> const& crQuat) {
	return crQuat.w;
}

template <Numeric T>
Quaternion<T> Mul(Quaternion<T> const& crLHS, Quaternion<T> const& crRHS) {
	ConstVectorView<T, 3> imLHS = Im(crLHS);
	ConstVectorView<T, 3> imRHS = Im(crRHS);

	T const& reLHS = Re(crLHS);
	T const& reRHS = Re(crRHS);

	return { reLHS * reRHS - Dot(imLHS, imRHS), Cross(imLHS, imRHS) + reLHS * imRHS + reRHS * imLHS };
}

template <Numeric T>
VectorView<T, 3> Im(Quaternion<T> const& crQuat) {

	// if the data of the quaternion is stored in contiguous memory...
	if constexpr (
		offsetof(Quaternion<T>, w) == 0 &&
		offsetof(Quaternion<T>, x) == sizeof(T) &&
		offsetof(Quaternion<T>, y) == 2 * sizeof(T) &&
		offsetof(Quaternion<T>, z) == 3 * sizeof(T)
		) {

		return VectorView<T, 3>(const_cast<T*>(&crQuat.x));
	}

	else { 
		throw std::runtime_error("Elements of Quaternion are NOT stored in"
			" contiguous memory");
	}
}


template <CanUnaryPlus VT>
Quaternion<unary_plus_t<VT>> operator+(Quaternion<VT> const& crQuat) {
    Quaternion<unary_plus_t<VT>> result;

    result.w = +crQuat.w;
    result.x = +crQuat.x;
    result.y = +crQuat.y;
    result.z = +crQuat.z;

    return result;
}

template <CanNegate VT>
Quaternion<negation_t<VT>> operator-(Quaternion<VT> const& crQuat) {
    Quaternion<negation_t<VT>> result;

    result.w = -crQuat.w;
    result.x = -crQuat.x;
    result.y = -crQuat.y;
    result.z = -crQuat.z;

    return result;
}

template <CanDereference VT>
Quaternion<dereference_t<VT>> operator*(Quaternion<VT> const& crQuat) {
    Quaternion<dereference_t<VT>> result;

    result.w = *crQuat.w;
    result.x = *crQuat.x;
    result.y = *crQuat.y;
    result.z = *crQuat.z;

    return result;
}

template <CanBitwiseNot VT>
Quaternion<bitwise_not_t<VT>> operator~(Quaternion<VT> const& crQuat) {

    Quaternion<bitwise_not_t<VT>> result;

    result.w = ~crQuat.w;
    result.x = ~crQuat.x;
    result.y = ~crQuat.y;
    result.z = ~crQuat.z;

    return result;
}

template <CanLogicalNot VT>
bool operator!(Quaternion<VT> const& crQuat) {
    if (!crQuat.w) {}
    else { return false; }

    if (!crQuat.x) {}
    else { return false; }

    if (!crQuat.y) {}
    else { return false; }

    if (!crQuat.z) {}
    else { return false; }

    return true;
}

template <CanPrefixIncrement VT>
Quaternion<VT>& operator++(Quaternion<VT>& rQuat) {
    ++rQuat.w;
    ++rQuat.x;
    ++rQuat.y;
    ++rQuat.z;

    return rQuat;
}

template <CanPostfixIncrement VT>
Quaternion<VT> operator++(Quaternion<VT>& rQuat, int) {
    Quaternion<VT> q = rQuat;

    rQuat.w++;
    rQuat.x++;
    rQuat.y++;
    rQuat.z++;

    return q;
}

template <CanPrefixDecrement VT>
Quaternion<VT>& operator--(Quaternion<VT>& rQuat) {
    --rQuat.w;
    --rQuat.x;
    --rQuat.y;
    --rQuat.z;

    return rQuat;
}

template <CanPostfixDecrement VT>
Quaternion<VT> operator--(Quaternion<VT>& rQuat, int) {
    Quaternion<VT> q = rQuat;

    rQuat.w--;
    rQuat.x--;
    rQuat.y--;
    rQuat.z--;

    return q;
}

template <typename VT1, typename VT2> requires CanSum<VT1, VT2>
Quaternion<sum_t<VT1, VT2>> operator+(Quaternion<VT1> const& crLHS, Quaternion<VT2> const& crRHS) {
    Quaternion<sum_t<VT1, VT2>> result;

    result.w = crLHS.w + crRHS.w;
    result.x = crLHS.x + crRHS.x;
    result.y = crLHS.y + crRHS.y;
    result.z = crLHS.z + crRHS.z;

    return result;
}

template <typename VT, typename ST> requires CanSum<VT, ST>
Quaternion<sum_t<VT, ST>> operator+(Quaternion<VT> const& crLHS, ST const& crRHS) {
    Quaternion<sum_t<VT, ST>> result;

    result.w = crLHS.w + crRHS;
    result.x = crLHS.x + crRHS;
    result.y = crLHS.y + crRHS;
    result.z = crLHS.z + crRHS;

    return result;
}

template <typename ST, typename VT> requires CanSum<ST, VT>
Quaternion<sum_t<ST, VT>> operator+(ST const& crLHS, Quaternion<VT> const& crRHS) {
    Quaternion<sum_t<ST, VT>> result;

    result.w = crLHS + crRHS.w;
    result.x = crLHS + crRHS.x;
    result.y = crLHS + crRHS.y;
    result.z = crLHS + crRHS.z;

    return result;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator+=(DB1<VT1, Ds...>&& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<sum_t<VT1, VT2>, VT1>, "Cannot apply compound-sum-assignment to two DataBlock when "
        "the sum type of the element types of the DataBlocks are not implicitly convertable to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] += crRHS[i];
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator+=(DB1<VT1, Ds...>& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<sum_t<VT1, VT2>, VT1>, "Cannot apply compound-sum-assignment to two DataBlock when "
        "the sum type of the element types of the DataBlocks are not implicitly convertable to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] += crRHS[i];
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator+=(DB<VT, Ds...>&& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<sum_t<VT, ST>, VT>, "Cannot apply compound-sum-assignment to a DataBlock and another type when "
        "the sum type of the right-hand type and the element type of the DataBlock is not implicitly convertable to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] += crRHS;
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator+=(DB<VT, Ds...>& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<sum_t<VT, ST>, VT>, "Cannot apply compound-sum-assignment to a DataBlock and another type when "
        "the sum type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] += crRHS;
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto operator-(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanSubtract<VT1, VT2>, "Cannot subtract DataBlocks whose element types cannot be subtracted.");

    DataBlock<difference_t<VT1, VT2>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] - crRHS[i];
    }

    return result;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
auto operator-(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanSubtract<VT, ST>, "Cannot subtract a DataBlock and some other type when"
        " the other type and the elements of the DataBlock cannot be subtracted.");

    DataBlock<difference_t<VT, ST>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] - crRHS;
    }

    return result;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
auto operator-(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanSum<ST, VT>, "Cannot subtract a DataBlock and some other type when"
        " the other type and the elements of the DataBlock cannot be subtracted.");

    DataBlock<difference_t<ST, VT>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS - crRHS[i];
    }

    return result;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator-=(DB1<VT1, Ds...>&& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<difference_t<VT1, VT2>, VT1>, "Cannot apply compound-subtraction-assignment to two DataBlock when "
        "the difference type of the element types of the DataBlocks are not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] -= crRHS[i];
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator-=(DB1<VT1, Ds...>& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<difference_t<VT1, VT2>, VT1>, "Cannot apply compound-subtraction-assignment to two DataBlocks when "
        "the difference type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] -= crRHS[i];
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator-=(DB<VT, Ds...>&& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<difference_t<VT, ST>, VT>, "Cannot apply compound-subtraction-assignment to a DataBlock and another type when "
        "the difference type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] -= crRHS;
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator-=(DB<VT, Ds...>& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<difference_t<VT, ST>, VT>, "Cannot apply compound-subtraction-assignment to a DataBlock and another type when "
        "the difference type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] -= crRHS;
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
DataBlock<product_t<VT1, VT2>, Ds...> operator*(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanMultiply<VT1, VT2>, "Cannot multiply DataBlocks whose element types cannot be multiplied.");

    DataBlock<product_t<VT1, VT2>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] * crRHS[i];
    }

    return result;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
auto operator*(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanMultiply<VT, ST>, "Cannot multiply a DataBlock and some other type when"
        " the other type and the elements of the DataBlock cannot be multiplied.");

    DataBlock<product_t<VT, ST>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] * crRHS;
    }

    return result;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
auto operator*(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanMultiply<ST, VT>, "Cannot multiply a DataBlock and some other type when"
        " the other type and the elements of the DataBlock cannot be multiplied.");

    DataBlock<product_t<ST, VT>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS * crRHS[i];
    }

    return result;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator*=(DB1<VT1, Ds...>&& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<product_t<VT1, VT2>, VT1>, "Cannot apply compound-multiplication-assignment to two DataBlocks when "
        "the product type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] *= crRHS[i];
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator*=(DB1<VT1, Ds...>& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<product_t<VT1, VT2>, VT1>, "Cannot apply compound-multiplication-assignment to two DataBlocks when "
        "the product type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] *= crRHS[i];
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator*=(DB<VT, Ds...>&& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<product_t<VT, ST>, VT>, "Cannot apply compound-multiplication-assignment to a DataBlock and another type when "
        "the product type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] *= crRHS;
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator*=(DB<VT, Ds...>& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<product_t<VT, ST>, VT>, "Cannot apply compound-multiplication-assignment to a DataBlock and another type when "
        "the product type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] *= crRHS;
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto operator/(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanDivide<VT1, VT2>, "Cannot divide DataBlocks whose element types cannot be divided.");

    DataBlock<quotient_t<VT1, VT2>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] / crRHS[i];
    }

    return result;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
auto operator/(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanDivide<VT, ST>, "Cannot divide a DataBlock and some other type when"
        " the other type and the elements of the DataBlock cannot be divided.");

    DataBlock<quotient_t<VT, ST>, Ds...> result;
    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] / crRHS;
    }

    return result;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
auto operator/(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanDivide<ST, VT>, "Cannot divide a DataBlock and some other type when"
        " the other type and the elements of the DataBlock cannot be divided.");

    DataBlock<quotient_t<ST, VT>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS / crRHS[i];
    }

    return result;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator/=(DB1<VT1, Ds...>&& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<quotient_t<VT1, VT2>, VT1>, "Cannot apply compound-division-assignment to two DataBlocks when "
        "the quotient type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] /= crRHS[i];
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator/=(DB1<VT1, Ds...>& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<quotient_t<VT1, VT2>, VT1>, "Cannot apply compound-division-assignment to two DataBlocks when "
        "the quotient type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] /= crRHS[i];
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator/=(DB<VT, Ds...>&& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<quotient_t<VT, ST>, VT>, "Cannot apply compound-division-assignment to a DataBlock and another type when "
        "the quotient type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] /= crRHS;
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator/=(DB<VT, Ds...>& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<quotient_t<VT, ST>, VT>, "Cannot apply compound-division-assignment to a DataBlock and another type when "
        "the quotient type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] /= crRHS;
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto operator%(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanModulus<VT1, VT2>, "Cannot take the modulus of DataBlocks for which the modulus of their element types cannot be taken.");

    DataBlock<modulus_t<VT1, VT2>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] % crRHS[i];
    }

    return result;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
auto operator%(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanModulus<VT, ST>, "Cannot take the modulus of a DataBlock and some other type when"
        " the modulus of the other type and the elements of the DataBlock cannot be taken.");

    DataBlock<modulus_t<VT, ST>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] % crRHS;
    }

    return result;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
auto operator%(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanModulus<ST, VT>, "Cannot take the modulus of a DataBlock and some other type when"
        " the modulus of the other type and the elements of the DataBlock cannot be taken.");

    DataBlock<modulus_t<ST, VT>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS % crRHS[i];
    }

    return result;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator%=(DB1<VT1, Ds...>&& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<modulus_t<VT1, VT2>, VT1>, "Cannot apply compound-modulus-assignment to two DataBlocks when "
        "the modulus type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] %= crRHS[i];
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator%=(DB1<VT1, Ds...>& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<modulus_t<VT1, VT2>, VT1>, "Cannot apply compound-modulus-assignment to two DataBlocks when "
        "the modulus type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] %= crRHS[i];
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator%=(DB<VT, Ds...>&& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<modulus_t<VT, ST>, VT>, "Cannot apply compound-modulus-assignment to a DataBlock and another type when "
        "the modulus type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] %= crRHS;
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator%=(DB<VT, Ds...>& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<modulus_t<VT, ST>, VT>, "Cannot apply compound-modulus-assignment to a DataBlock and another type when "
        "the modulus type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] %= crRHS;
    }

    return rLHS;
}


template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto operator^(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanBitwiseXor<VT1, VT2>, "Cannot apply bitwise xor to DataBlocks for which bitwise xor cannot be applied to their element types.");

    DataBlock<bitwise_xor_t<VT1, VT2>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] ^ crRHS[i];
    }

    return result;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
auto operator^(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanBitwiseXor<VT, ST>, "Cannot apply bitwise xor to a DataBlock and some other type when"
        " bitwise xor cannot be applied to the other type and the elements of the DataBlock.");

    DataBlock<bitwise_xor_t<VT, ST>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] ^ crRHS;
    }

    return result;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
auto operator^(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanBitwiseXor<ST, VT>, "Cannot apply bitwise xor to a DataBlock and some other type when"
        " bitwise xor cannot be applied to the other type and the elements of the DataBlock.");

    DataBlock<bitwise_xor_t<ST, VT>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS ^ crRHS[i];
    }

    return result;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator^=(DB1<VT1, Ds...>&& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<bitwise_xor_t<VT1, VT2>, VT1>, "Cannot apply compound-bitwise-xor-assignment to two DataBlocks when "
        "the bitwise-xor type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] ^= crRHS[i];
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator^=(DB1<VT1, Ds...>& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<bitwise_xor_t<VT1, VT2>, VT1>, "Cannot apply compound-bitwise-xor-assignment to two DataBlocks when "
        "the bitwise-xor type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] ^= crRHS[i];
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator^=(DB<VT, Ds...>&& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<bitwise_xor_t<VT, ST>, VT>, "Cannot apply compound-bitwise-xor-assignment to a DataBlock and another type when "
        "the bitwise-xor type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] ^= crRHS;
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator^=(DB<VT, Ds...>& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<bitwise_xor_t<VT, ST>, VT>, "Cannot apply compound-bitwise-xor-assignment to a DataBlock and another type when "
        "the bitwise-xor type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] ^= crRHS;
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto operator&(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanBitwiseAnd<VT1, VT2>, "Cannot apply bitwise and to DataBlocks for which bitwise and cannot be applied to their element types.");

    DataBlock<bitwise_and_t<VT1, VT2>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] & crRHS[i];
    }

    return result;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
auto operator&(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanBitwiseAnd<VT, ST>, "Cannot apply bitwise and to a DataBlock and some other type when"
        " bitwise and cannot be applied to the other type and the elements of the DataBlock.");

    DataBlock<bitwise_and_t<VT, ST>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] & crRHS;
    }

    return result;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
auto operator&(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanBitwiseAnd<ST, VT>, "Cannot apply bitwise and to a DataBlock and some other type when"
        " bitwise and cannot be applied to the other type and the elements of the DataBlock.");

    DataBlock<bitwise_and_t<ST, VT>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS & crRHS[i];
    }

    return result;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator&=(DB1<VT1, Ds...>&& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<bitwise_and_t<VT1, VT2>, VT1>, "Cannot apply compound-bitwise-and-assignment to two DataBlocks when "
        "the bitwise-and type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] &= crRHS[i];
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator&=(DB1<VT1, Ds...>& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<bitwise_and_t<VT1, VT2>, VT1>, "Cannot apply compound-bitwise-and-assignment to two DataBlocks when "
        "the bitwise-and type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] &= crRHS[i];
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator&=(DB<VT, Ds...>&& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<bitwise_and_t<VT, ST>, VT>, "Cannot apply compound-bitwise-and-assignment to a DataBlock and another type when "
        "the bitwise-and type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] &= crRHS;
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator&=(DB<VT, Ds...>& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<bitwise_and_t<VT, ST>, VT>, "Cannot apply compound-bitwise-and-assignment to a DataBlock and another type when "
        "the bitwise-and type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] &= crRHS;
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto operator|(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanBitwiseOr<VT1, VT2>, "Cannot apply bitwise-or to DataBlocks for which bitwise-or cannot be applied to their element types.");

    DataBlock<bitwise_or_t<VT1, VT2>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] | crRHS[i];
    }

    return result;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
auto operator|(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanBitwiseOr<VT, ST>, "Cannot apply bitwise-or to a DataBlock and some other type when"
        " bitwise-or cannot be applied to the other type and the elements of the DataBlock.");

    DataBlock<bitwise_or_t<VT, ST>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] | crRHS;
    }

    return result;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
auto operator|(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanBitwiseOr<ST, VT>, "Cannot apply bitwise-or to a DataBlock and some other type when"
        " bitwise-or cannot be applied to the other type and the elements of the DataBlock.");

    DataBlock<bitwise_or_t<ST, VT>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS | crRHS[i];
    }

    return result;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator|=(DB1<VT1, Ds...>&& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<bitwise_or_t<VT1, VT2>, VT1>, "Cannot apply compound-bitwise-or-assignment to two DataBlocks when "
        "the bitwise-or type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] |= crRHS[i];
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator|=(DB1<VT1, Ds...>& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<bitwise_or_t<VT1, VT2>, VT1>, "Cannot apply compound-bitwise-or-assignment to two DataBlocks when "
        "the bitwise-or type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] |= crRHS[i];
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator|=(DB<VT, Ds...>&& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<bitwise_or_t<VT, ST>, VT>, "Cannot apply compound-bitwise-or-assignment to a DataBlock and another type when "
        "the bitwise-or type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] |= crRHS;
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator|=(DB<VT, Ds...>& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<bitwise_or_t<VT, ST>, VT>, "Cannot apply compound-bitwise-or-assignment to a DataBlock and another type when "
        "the bitwise-or type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] |= crRHS;
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
bool operator<(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanLessThan<VT1, VT2>, "Cannot apply a less-than comparison to DataBlocks whose element types cannot be tested in less-than comparisons.");
    static_assert(IsExplicitlyConvertible<less_than_t<VT1, VT2>, bool>, "Cannot apply a less-than comparison to DataBlocks whose element types yield a type that "
        "cannot be explicitly converted to bool under less-than comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] < crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
bool operator<(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanLessThan<VT, ST>, "Cannot apply a less-than comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in less-than comparisons.");
    static_assert(IsExplicitlyConvertible<less_than_t<VT, ST>, bool>, "Cannot apply a less-than comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under less-than comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] < crRHS) {}
        else { return false; }
    }

    return true;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
bool operator<(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanLessThan<ST, VT>, "Cannot apply a less-than comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in less-than comparisons.");
    static_assert(IsExplicitlyConvertible<less_than_t<ST, VT>, bool>, "Cannot apply a less-than comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under less-than comparison.");

    for (int i = 0; i < crRHS.nelems; ++i) {
        if (crLHS < crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
bool operator>(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanGreaterThan<VT1, VT2>, "Cannot apply a greater-than comparison to DataBlocks whose element types cannot be tested in greater-than comparisons.");
    static_assert(IsExplicitlyConvertible<greater_than_t<VT1, VT2>, bool>, "Cannot apply a greater-than comparison to DataBlocks whose element types yield a type that "
        "cannot be explicitly converted to bool under greater-than comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] > crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
bool operator>(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanGreaterThan<VT, ST>, "Cannot apply a greater-than comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in greater-than comparisons.");
    static_assert(IsExplicitlyConvertible<greater_than_t<VT, ST>, bool>, "Cannot apply a greater-than comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under greater-than comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] > crRHS) {}
        else { return false; }
    }

    return true;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
bool operator>(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanGreaterThan<ST, VT>, "Cannot apply a greater-than comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in greater-than comparisons.");
    static_assert(IsExplicitlyConvertible<greater_than_t<ST, VT>, bool>, "Cannot apply a greater-than comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under greater-than comparison.");

    for (int i = 0; i < crRHS.nelems; ++i) {
        if (crLHS > crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto operator<<(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanLeftShift<VT1, VT2>, "Cannot left-shift DataBlocks whose element types cannot be left-shifted.");

    DataBlock<left_shift_t<VT1, VT2>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] << crRHS[i];
    }

    return result;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
auto operator<<(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanLeftShift<VT, ST>, "Cannot left-shift a DataBlock and some other type when"
        " the other type and the elements of the DataBlock cannot be left-shifted.");

    DataBlock<left_shift_t<VT, ST>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] << crRHS;
    }

    return result;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
auto operator<<(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanLeftShift<ST, VT>, "Cannot left-shift a DataBlock and some other type when"
        " the other type and the elements of the DataBlock cannot be left-shifted.");

    DataBlock<left_shift_t<ST, VT>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS << crRHS[i];
    }

    return result;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator<<=(DB1<VT1, Ds...>&& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<left_shift_t<VT1, VT2>, VT1>, "Cannot apply compound-left-shift-assignment to two DataBlocks when "
        "the left-shift type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] <<= crRHS[i];
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator<<=(DB1<VT1, Ds...>& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<left_shift_t<VT1, VT2>, VT1>, "Cannot apply compound-left-shift-assignment to two DataBlocks when "
        "the left-shift type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] <<= crRHS[i];
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator<<=(DB<VT, Ds...>&& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<left_shift_t<VT, ST>, VT>, "Cannot apply compound-left-shift-assignment to a DataBlock and another type when "
        "the left-shift type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] <<= crRHS;
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator<<=(DB<VT, Ds...>& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<left_shift_t<VT, ST>, VT>, "Cannot apply compound-left-shift-assignment to a DataBlock and another type when "
        "the left-shift type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] <<= crRHS;
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto operator>>(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanRightShift<VT1, VT2>, "Cannot right-shift DataBlocks whose element types cannot be right-shifted.");

    DataBlock<right_shift_t<VT1, VT2>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] >> crRHS[i];
    }

    return result;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
auto operator>>(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanRightShift<VT, ST>, "Cannot right-shift a DataBlock and some other type when"
        " the other type and the elements of the DataBlock cannot be right-shifted.");

    DataBlock<right_shift_t<VT, ST>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] >> crRHS;
    }

    return result;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
auto operator>>(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanRightShift<ST, VT>, "Cannot right-shift a DataBlock and some other type when"
        " the other type and the elements of the DataBlock cannot be right-shifted.");

    DataBlock<right_shift_t<ST, VT>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS >> crRHS[i];
    }

    return result;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator>>=(DB1<VT1, Ds...>&& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<right_shift_t<VT1, VT2>, VT1>, "Cannot apply compound-right-shift-assignment to two DataBlocks when "
        "the right-shift type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] >>= crRHS[i];
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator>>=(DB1<VT1, Ds...>& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<right_shift_t<VT1, VT2>, VT1>, "Cannot apply compound-right-shift-assignment to two DataBlocks when "
        "the right-shift type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] >>= crRHS[i];
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator>>=(DB<VT, Ds...>&& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<right_shift_t<VT, ST>, VT>, "Cannot apply compound-right-shift-assignment to a DataBlock and another type when "
        "the right-shift type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] >>= crRHS;
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator>>=(DB<VT, Ds...>& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<right_shift_t<VT, ST>, VT>, "Cannot apply compound-right-shift-assignment to a DataBlock and another type when "
        "the right-shift type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] >>= crRHS;
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
bool operator==(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanLogicalEquality<VT1, VT2>, "Cannot apply a logical-equality comparison to DataBlocks whose element types cannot be tested in logical-equality comparisons.");
    static_assert(IsExplicitlyConvertible<logical_equality_t<VT1, VT2>, bool>, "Cannot apply a logical-equality comparison to DataBlocks whose element types yield a type that "
        "cannot be explicitly converted to bool under logical-equality comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] == crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
bool operator==(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanLogicalEquality<VT, ST>, "Cannot apply a logical-equality comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in logical-equality comparisons.");
    static_assert(IsExplicitlyConvertible<logical_equality_t<VT, ST>, bool>, "Cannot apply a logical-equality comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under logical-equality comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] == crRHS) {}
        else { return false; }
    }

    return true;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
bool operator==(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanLogicalEquality<ST, VT>, "Cannot apply a logical-equality comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in logical-equality comparisons.");
    static_assert(IsExplicitlyConvertible<logical_equality_t<ST, VT>, bool>, "Cannot apply a logical-equality comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under logical-equality comparison.");

    for (int i = 0; i < crRHS.nelems; ++i) {
        if (crLHS == crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
bool operator!=(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanLogicalInequality<VT1, VT2>, "Cannot apply a logical-inequality comparison to DataBlocks whose element types cannot be tested in logical-inequality comparisons.");
    static_assert(IsExplicitlyConvertible<logical_inequality_t<VT1, VT2>, bool>, "Cannot apply a logical-inequality comparison to DataBlocks whose element types yield a type that "
        "cannot be explicitly converted to bool under logical-inequality comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] != crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
bool operator!=(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanLogicalInequality<VT, ST>, "Cannot apply a logical-inequality comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in logical-inequality comparisons.");
    static_assert(IsExplicitlyConvertible<logical_inequality_t<VT, ST>, bool>, "Cannot apply a logical-inequality comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under logical-inequality comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] != crRHS) {}
        else { return false; }
    }

    return true;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
bool operator!=(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanLogicalInequality<ST, VT>, "Cannot apply a logical-inequality comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in logical-inequality comparisons.");
    static_assert(IsExplicitlyConvertible<logical_inequality_t<ST, VT>, bool>, "Cannot apply a logical-inequality comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under logical-inequality comparison.");

    for (int i = 0; i < crRHS.nelems; ++i) {
        if (crLHS != crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
bool operator<=(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanLessThanOrEqual<VT1, VT2>, "Cannot apply a less-than-or-equal comparison to DataBlocks whose element types cannot be tested in less-than-or-equal comparisons.");
    static_assert(IsExplicitlyConvertible<less_than_or_equal_t<VT1, VT2>, bool>, "Cannot apply a less-than-or-equal comparison to DataBlocks whose element types yield a type that "
        "cannot be explicitly converted to bool under less-than-or-equal comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] <= crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
bool operator<=(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanLessThanOrEqual<VT, ST>, "Cannot apply a less-than-or-equal comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in less-than-or-equal comparisons.");
    static_assert(IsExplicitlyConvertible<less_than_or_equal_t<VT, ST>, bool>, "Cannot apply a less-than-or-equal comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under less-than-or-equal comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] <= crRHS) {}
        else { return false; }
    }

    return true;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
bool operator<=(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanLessThanOrEqual<ST, VT>, "Cannot apply a less-than-or-equal comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in less-than-or-equal comparisons.");
    static_assert(IsExplicitlyConvertible<less_than_or_equal_t<ST, VT>, bool>, "Cannot apply a less-than-or-equal comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under less-than-or-equal comparison.");

    for (int i = 0; i < crRHS.nelems; ++i) {
        if (crLHS <= crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
bool operator>=(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanGreaterThanOrEqual<VT1, VT2>, "Cannot apply a greater-than-or-equal comparison to DataBlocks whose element types cannot be tested in greater-than-or-equal comparisons.");
    static_assert(IsExplicitlyConvertible<greater_than_or_equal_t<VT1, VT2>, bool>, "Cannot apply a greater-than-or-equal comparison to DataBlocks whose element types yield a type that "
        "cannot be explicitly converted to bool under greater-than-or-equal comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] >= crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
bool operator>=(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanGreaterThanOrEqual<VT, ST>, "Cannot apply a greater-than-or-equal comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in greater-than-or-equal comparisons.");
    static_assert(IsExplicitlyConvertible<greater_than_or_equal_t<VT, ST>, bool>, "Cannot apply a greater-than-or-equal comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under greater-than-or-equal comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] >= crRHS) {}
        else { return false; }
    }

    return true;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
bool operator>=(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanGreaterThanOrEqual<ST, VT>, "Cannot apply a greater-than-or-equal comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in greater-than-or-equal comparisons.");
    static_assert(IsExplicitlyConvertible<greater_than_or_equal_t<ST, VT>, bool>, "Cannot apply a greater-than-or-equal comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under greater-than-or-equal comparison.");

    for (int i = 0; i < crRHS.nelems; ++i) {
        if (crLHS >= crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
three_way_compare_t<VT1, VT2> operator<=>(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanThreeWayCompare<VT1, VT2>, "Cannot apply a three-way comparison to DataBlocks whose element types cannot be tested in three-way comparisons.");
    static_assert(IsValidThreeWayCompareReturn<three_way_compare_t<VT1, VT2>>, "Cannot apply a three-way comparison to DataBlocks whose element types yield a type that "
        "is not a canonical three-way comparison return type.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        three_way_compare_t<VT1, VT2> result = crLHS[i] <=> crRHS[i];
        if (result > 0) return result;
        if (result < 0) return result;
    }

    return crLHS[crLHS.nelems - 1] <=> crRHS[crRHS.nelems - 1];
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
three_way_compare_t<VT, ST> operator<=>(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanThreeWayCompare<VT, ST>, "Cannot apply a three-way comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in three-way comparisons.");
    static_assert(IsValidThreeWayCompareReturn<three_way_compare_t<VT, ST>>, "Cannot apply a three-way comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that is not a canonical "
        "three-way comparison return type.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        three_way_compare_t<VT, ST> result = crLHS[i] <=> crRHS;
        if (result > 0) return result;
        if (result < 0) return result;
    }

    return crLHS[crLHS.nelems - 1] <=> crRHS;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
three_way_compare_t<ST, VT> operator<=>(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanThreeWayCompare<ST, VT>, "Cannot apply a three-way comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in three-way comparisons.");
    static_assert(IsValidThreeWayCompareReturn<three_way_compare_t<ST, VT>>, "Cannot apply a three-way comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that is not a canonical three-way"
        " comparison return type.");

    for (int i = 0; i < crRHS.nelems; ++i) {
        three_way_compare_t<ST, VT> result = crLHS <=> crRHS[i];
        if (result > 0) return result;
        if (result < 0) return result;
    }

    return crLHS <=> crRHS[crRHS.nelems - 1];
}

 quaternion operations w/ Vectors (i.e., single-dimensional DataBlocks)
template <Numeric T>
Vector<T, 4> QuatMul(Vector<T, 4> const& crLHS, Vector<T, 4> const& crRHS) {
	ConstVectorView<T, 3> imLHS = QuatIm(crLHS);
	ConstVectorView<T, 3> imRHS = QuatIm(crRHS);

	T const& wLHS = QuatRe(crLHS);
	T const& wRHS = QuatRe(crRHS);

	return { wLHS * wRHS - Dot(imLHS, imRHS), Cross(imLHS, imRHS) + wLHS * imRHS + wRHS * imLHS };
}

template <Numeric T>
T QuatSquaredMagnitude(Vector<T, 4> const& crQuat) {
	return SquaredMagnitude(crQuat);
}

template <Numeric T>
T QuatMagnitude(Vector<T, 4> const& crQuat) {
	return Magnitude(crQuat);
}

template <Numeric T>
T& QuatRe(Vector<T, 4>& rQuat) {
	return rQuat[0];
}

template <Numeric T>
T const& QuatRe(Vector<T, 4> const& crQuat) {
	return crQuat[0];
}

template <Numeric T>
VectorView<T, 3> QuatIm(Vector<T, 4>& rQuat) {
	return rQuat[1_idx >> 4_idx];
}

template <Numeric T>
ConstVectorView<T, 3> QuatIm(Vector<T, 4> const& crQuat) {
	return crQuat[1_idx >> 4_idx];
}

template <Numeric T>
Vector<T, 4> QuatConj(Vector<T, 4> const& crVec) {
	return {crVec[0], -crVec[1], -crVec[2], -crVec[3]};
}

template <Numeric T>
Vector<T, 4> QuatInv(Vector<T, 4> const& crQuat) {
	return QuatConj(crQuat) / QuatSquaredMagnitude(crQuat);
}

using quat  = Quaternion<float>;
using iquat = Quaternion<int>;
using uquat = Quaternion<unsigned int>;
using dquat = Quaternion<double>;

#endif