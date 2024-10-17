#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <concepts>
#include <numbers>
#include <string>

#include "data-block.hpp"
#include "vector.hpp"

enum Axis {
	X = 0, Y = 1, Z = 2
};

template <typename T, size_t R, size_t C>
using Matrix = DataBlock<T, R, C>;

template <size_t R, size_t C>
using mat = Matrix<float, R, C>;

using mat2 = Matrix<float, 2, 2>;
using mat2x2 = Matrix<float, 2, 2>;
using mat2x3 = Matrix<float, 2, 3>;
using mat2x4 = Matrix<float, 2, 4>;

using mat3x2 = Matrix<float, 3, 2>;
using mat3 = Matrix<float, 3, 3>;
using mat3x3 = Matrix<float, 3, 3>;
using mat3x4 = Matrix<float, 3, 4>;

using mat4x2 = Matrix<float, 4, 2>;
using mat4x3 = Matrix<float, 4, 3>;
using mat4 = Matrix<float, 4, 4>;
using mat4x4 = Matrix<float, 4, 4>;

template <size_t R, size_t C>
using dmat = Matrix<double, R, C>;

using dmat2 = Matrix<double, 2, 2>;
using dmat2x2 = Matrix<double, 2, 2>;
using dmat2x3 = Matrix<double, 2, 3>;
using dmat2x4 = Matrix<double, 2, 4>;

using dmat3x2 = Matrix<double, 3, 2>;
using dmat3 = Matrix<double, 3, 3>;
using dmat3x3 = Matrix<double, 3, 3>;
using dmat3x4 = Matrix<double, 3, 4>;

using dmat4x2 = Matrix<double, 4, 2>;
using dmat4x3 = Matrix<double, 4, 3>;
using dmat4 = Matrix<double, 4, 4>;
using dmat4x4 = Matrix<double, 4, 4>;

template <typename T, size_t R, size_t C>
using MatrixView = DataBlockView<T, R, C>;

template <size_t R, size_t C>
using mat_view = MatrixView<float, R, C>;

using mat2_view = MatrixView<float, 2, 2>;
using mat2x2_view = MatrixView<float, 2, 2>;
using mat2x3_view = MatrixView<float, 2, 3>;
using mat2x4_view = MatrixView<float, 2, 4>;

using mat3x2_view = MatrixView<float, 3, 2>;
using mat3_view = MatrixView<float, 3, 3>;
using mat3x3_view = MatrixView<float, 3, 3>;
using mat3x4_view = MatrixView<float, 3, 4>;

using mat4x2_view = MatrixView<float, 4, 2>;
using mat4x3_view = MatrixView<float, 4, 3>;
using mat4_view = MatrixView<float, 4, 4>;
using mat4x4_view = MatrixView<float, 4, 4>;

template <size_t R, size_t C>
using dmat_view = MatrixView<double, R, C>;

using dmat2_view = MatrixView<double, 2, 2>;
using dmat2x2_view = MatrixView<double, 2, 2>;
using dmat2x3_view = MatrixView<double, 2, 3>;
using dmat2x4_view = MatrixView<double, 2, 4>;

using dmat3x2_view = MatrixView<double, 3, 2>;
using dmat3_view = MatrixView<double, 3, 3>;
using dmat3x3_view = MatrixView<double, 3, 3>;
using dmat3x4_view = MatrixView<double, 3, 4>;

using dmat4x2_view = MatrixView<double, 4, 2>;
using dmat4x3_view = MatrixView<double, 4, 3>;
using dmat4_view = MatrixView<double, 4, 4>;
using dmat4x4_view = MatrixView<double, 4, 4>;

template <Numeric T, size_t R1, size_t C1R2, size_t C2>
Matrix<T, R1, C2> Mul(Matrix<T, R1, C1R2> const& crLHS, Matrix<T, C1R2, C2> const& crRHS) {
	Matrix<T, R1, C2> result({.fill = 0});

	for (int i = 0; i < R1; ++i) {
		for (int j = 0; j < C2; ++j) {
			for (int k = 0; k < C1R2; ++k) {
				result(i, j) += crLHS(i, k) * crRHS(k, j);
			}
		}
	}

	return result;
}

template <Numeric T, size_t R, size_t C>
Vector<T, R> Mul(Matrix<T, R, C> const& crLHS, Vector<T, C> const& crRHS) {
	Vector<T, R> result({.fill = 0});

	for (int i = 0; i < R; ++i) {
		for (int j = 0; j < C; ++j) {
			result(i) += crLHS(i, j) * crRHS(j);
		}
	}

	return result;
}

template <Numeric T, size_t R, size_t C>
Vector<T, C> Mul(Vector<T, R> const& crLHS, Matrix<T, R, C> const& crRHS) {
	Vector<T, C> result({.fill = 0});

	for (int i = 0; i < R; ++i) {
		for (int j = 0; j < C; ++j) {
			result(i) += crLHS(i) * crRHS(i, j);
		}
	}

	return result;
}

template <Numeric T>
Matrix<T, 4, 4> TranslationMatrix(T x, T y, T z) {
	Matrix<T, 4, 4> result = Matrix<T, 4, 4>::Identity();

	result[3]  = x; // result(0, 3)
	result[7]  = y; // result(1, 3)
	result[11] = z; // result(2, 3)

	return result;
}

template <Numeric T, size_t R = 4, size_t C = 4>
requires ((R == 4 && C == 4) || (R == 2 && C == 2))
Matrix<T, R, C> XRotationMatrix(T radians) {
	if constexpr (R == 4 && C == 4) {
		return {
			1, 0,             0,            0,
			0, cos(radians), -sin(radians), 0,
			0, sin(radians),  cos(radians), 0,
			0, 0,             0,            1
		};
	}

	else { // R == 2 && C == 2
		return {
			cos(radians), -sin(radians),
			sin(radians),  cos(radians)
		};
	}
}

template <Numeric T, size_t R = 4, size_t C = 4>
requires ((R == 4 && C == 4) || (R == 2 && C == 2))
Matrix<T, R, C> YRotationMatrix(T radians) {
	if constexpr (R == 4 && C == 4) {
		return {
			 cos(radians), 0, sin(radians), 0,
		     0,            1, 0,            0,
			-sin(radians), 0, cos(radians), 0,
			 0,            0, 0,            1
		};
	}

	else { // R == 2 && C == 2
		return {
			 cos(radians),  sin(radians),
			-sin(radians),  cos(radians)
		};
	}
}

template <Numeric T, size_t R = 4, size_t C = 4>
requires ((R == 4 && C == 4) || (R == 2 && C == 2))
Matrix<T, R, C> ZRotationMatrix(T radians) {
	if constexpr (R == 4 && C == 4) {
		return {
			 cos(radians), -sin(radians), 0, 0,
			 sin(radians),  cos(radians), 0, 0,
			 0,             0,            1, 0,
			 0,             0,            0, 1
		};
	}

	else { // R == 2 && C == 2
		return {
			 cos(radians),  -sin(radians),
			 sin(radians),   cos(radians)
		};
	}
}

template <Numeric T>
Matrix<T, 4, 4> ScaleMatrix(T x, T y, T z) {
	return {
		x, 0, 0, 0,
		0, y, 0, 0,
		0, 0, z, 0,
		0, 0, 0, 1
	};
}

template <Numeric T>
void SetTranslation(T x, T y, T z, Matrix<T, 4, 4>& rMat) {
	rMat[3]  = x; // rMat(0, 3)
	rMat[7]  = y; // rMat(1, 3)
	rMat[11] = z; // rMat(2, 3)
}

template <Numeric T>
void SetTranslation(Vector<T, 3> const& crVec, Matrix<T, 4, 4>& rMat) {
	rMat[3]  = crVec[0]; // rMat(0, 3)
	rMat[7]  = crVec[1]; // rMat(1, 3)
	rMat[11] = crVec[2]; // rMat(2, 3)
}

template <Numeric T>
void SetTranslation(VectorView<T, 3> const& crVecView, Matrix<T, 4, 4>& rMat) {
	rMat[3]  = crVecView[0]; // rMat(0, 3)
	rMat[7]  = crVecView[1]; // rMat(1, 3)
	rMat[11] = crVecView[2]; // rMat(2, 3)
}

template <Numeric T>
void SetXTranslation(T x, Matrix<T, 4, 4>& rMat) {
	rMat[3]  = x; // rMat(0, 3)
}

template <Numeric T>
void SetYTranslation(T y, Matrix<T, 4, 4>& rMat) {
	rMat[7]  = y; // rMat(1, 3)
}

template <Numeric T>
void SetZTranslation(T z, Matrix<T, 4, 4>& rMat) {
	rMat[11] = z; // rMat(2, 3)
}

template <Numeric T>
void LTranslate(T x, T y, T z, Matrix<T, 4, 4>& rMat) {
	rMat[3]  += x; // rMat(0, 3)
	rMat[7]  += y; // rMat(1, 3)
	rMat[11] += z; // rMat(2, 3)
}

template <Numeric T>
void LXTranslate(T x, Matrix<T, 4, 4>& rMat) {
	rMat[3]  += x; // rMat(0, 3)
}

template <Numeric T>
void LYTranslate(T y, Matrix<T, 4, 4>& rMat) {
	rMat[7]  += y; // rMat(1, 3)
}

template <Numeric T>
void LZTranslate(T z, Matrix<T, 4, 4>& rMat) {
	rMat[11] += z; // rMat(2, 3)
}

template <Numeric T>
void RTranslate(T x, T y, T z, Matrix<T, 4, 4>& rMat) {
	Matrix<T, 4, 4> translator = TranslationMatrix(x, y, z);
	rMat = Mul(rMat, translator);
}

template <Numeric T>
void RXTranslate(T x, Matrix<T, 4, 4>& rMat) {
	Matrix<T, 4, 4> translator = TranslationMatrix(x, 0, 0);
	rMat = Mul(rMat, translator);
}

template <Numeric T>
void RYTranslate(T y, Matrix<T, 4, 4>& rMat) {
	Matrix<T, 4, 4> translator = TranslationMatrix(0, y, 0);
	rMat = Mul(rMat, translator);
}

template <Numeric T>
void RZTranslate(T z, Matrix<T, 4, 4>& rMat) {
	Matrix<T, 4, 4> translator = TranslationMatrix(0, 0, z);
	rMat = Mul(rMat, translator);
}

template <Numeric T>
void SetXRotation(T radians, Matrix<T, 4, 4>& rMat) {
	rMat[5]  =  cos(radians); // rMat(1, 1)
	rMat[9]  =  sin(radians); // rMat(2, 1)
	rMat[6]  = -sin(radians); // rMat(1, 2)
	rMat[10] =  cos(radians); // rMat(2, 2)
}

template <Numeric T>
void LXRotate(T radians, Matrix<T, 4, 4>& rMat) {
	Matrix<T, 4, 4> rotator = XRotationMatrix(radians);
	rMat = Mul(rotator, rMat);
}

template <Numeric T>
void RXRotate(T radians, Matrix<T, 4, 4>& rMat) {
	Matrix<T, 4, 4> rotator = XRotationMatrix(radians);
	rMat = Mul(rMat, rotator);
}

template <Numeric T>
void SetYRotation(T radians, Matrix<T, 4, 4>& rMat) {
	rMat[0]  =  cos(radians); // rMat(0, 0)
	rMat[8]  = -sin(radians); // rMat(2, 0)
	rMat[2]  =  sin(radians); // rMat(0, 2)
	rMat[10] =  cos(radians); // rMat(2, 2)
}

template <Numeric T>
void LYRotate(T radians, Matrix<T, 4, 4>& rMat) {
	Matrix<T, 4, 4> rotator = YRotationMatrix(radians);
	rMat = Mul(rotator, rMat);
}

template <Numeric T>
void RYRotate(T radians, Matrix<T, 4, 4>& rMat) {
	Matrix<T, 4, 4> rotator = YRotationMatrix(radians);
	rMat = Mul(rMat, rotator);
}

template <Numeric T>
void SetZRotation(T radians, Matrix<T, 4, 4>& rMat) {
	rMat[0] =  cos(radians); // rMat(0, 0)
	rMat[4] =  sin(radians); // rMat(1, 0)
	rMat[1] = -sin(radians); // rMat(0, 1)
	rMat[5] =  cos(radians); // rMat(1, 1)
}

template <Numeric T>
void LZRotate(T radians, Matrix<T, 4, 4>& rMat) {
	Matrix<T, 4, 4> rotator = ZRotationMatrix(radians);
	rMat = Mul(rotator, rMat);
}

template <Numeric T>
void RZRotate(T radians, Matrix<T, 4, 4>& rMat) {
	Matrix<T, 4, 4> rotator = ZRotationMatrix(radians);
	rMat = Mul(rMat, rotator);
}

template <Numeric T>
void SetScale(T x, T y, T z, Matrix<T, 4, 4>& rMat) {
	rMat[0]  = x; // rMat(0, 0)
	rMat[5]  = y; // rMat(1, 1)
	rMat[10] = z; // rMat(2, 2)
}

template <Numeric T>
void SetXScale(T x, Matrix<T, 4, 4>& rMat) {
	rMat[0]  = x; // rMat(0, 0)
}

template <Numeric T>
void SetYScale(T y, Matrix<T, 4, 4>& rMat) {
	rMat[5]  = y; // rMat(1, 1)
}

template <Numeric T>
void SetZScale(T z, Matrix<T, 4, 4>& rMat) {
	rMat[10] = z; // rMat(2, 2)
}

template <Numeric T>
void LScale(T x, T y, T z, Matrix<T, 4, 4>& rMat) {
	Matrix<T, 4, 4> scaler = ScaleMatrix(x, y, z);
	rMat = Mul(scaler, rMat);
}

template <Numeric T>
void LXScale(T x, Matrix<T, 4, 4>& rMat) {
	Matrix<T, 4, 4> scaler = ScaleMatrix(x, 1, 1);
	rMat = Mul(scaler, rMat);
}

template <Numeric T>
void LYScale(T y, Matrix<T, 4, 4>& rMat) {
	Matrix<T, 4, 4> scaler = ScaleMatrix(1, y, 1);
	rMat = Mul(scaler, rMat);
}

template <Numeric T>
void LZScale(T z, Matrix<T, 4, 4>& rMat) {
	Matrix<T, 4, 4> scaler = ScaleMatrix(1, 1, z);
	rMat = Mul(scaler, rMat);
}

template <Numeric T>
void RScale(T x, T y, T z, Matrix<T, 4, 4>& rMat) {
	Matrix<T, 4, 4> scaler = ScaleMatrix(x, y, z);
	rMat = Mul(rMat, scaler);
}

template <Numeric T>
void RXScale(T x, Matrix<T, 4, 4>& rMat) {
	Matrix<T, 4, 4> scaler = ScaleMatrix(x, 1, 1);
	rMat = Mul(rMat, scaler);
}

template <Numeric T>
void RYScale(T y, Matrix<T, 4, 4>& rMat) {
	Matrix<T, 4, 4> scaler = ScaleMatrix(1, y, 1);
	rMat = Mul(rMat, scaler);
}

template <Numeric T>
void RZScale(T z, Matrix<T, 4, 4>& rMat) {
	Matrix<T, 4, 4> scaler = ScaleMatrix(1, 1, z);
	rMat = Mul(rMat, scaler);
}


template <size_t R = 4, size_t C = 4> requires ((R == 4 && C == 4) || (R == 3 && C == 3))
auto NonUnitQuat2RotationMatrix(VectorType<4> auto const& crQuat) -> Matrix<ValueType<decltype(crQuat)>, R, C> {
	using T = ValueType<decltype(crQuat)>;

	T s = static_cast<T>(2) / QuatSquaredMagnitude(crQuat);
	T w = crQuat[0];
	T x = crQuat[1];
	T y = crQuat[2];
	T z = crQuat[3];

	if constexpr (R == 4 && C == 4) {
		return {
			1 - s * (y * y + z * z), s * (x * y - w * z),     s * (x * z + w * y),     0,
			s * (x * y + w * z),     1 - s * (x * x + z * z), s * (y * z - w * x),     0,
			s * (x * z - w * y),     s * (y * z + w * x),     1 - s * (x * x + y * y), 0,
			0,                       0,                       0,                       1
		};
	}

	else { // R == 3 && C == 3
		return {
			1 - s * (y * y + z * z), s * (x * y - w * z),     s * (x * z + w * y),
			s * (x * y + w * z),     1 - s * (x * x + z * z), s * (y * z - w * x),
			s * (x * z - w * y),     s * (y * z + w * x),     1 - s * (x * x + y * y)
		};
	}
}

template <size_t R = 4, size_t C = 4> requires ((R == 4 && C == 4) || (R == 3 && C == 3))
auto UnitQuat2RotationMatrix(VectorType<4> auto const& crQuat) -> Matrix<ValueType<decltype(crQuat)>, R, C> {
	using T = ValueType<decltype(crQuat)>;

	T w = crQuat[0];
	T x = crQuat[1];
	T y = crQuat[2];
	T z = crQuat[3];

	if constexpr (R == 4 && C == 4) {
		return {
			1 - 2 * (y * y + z * z), 2 * (x * y - w * z),     2 * (x * z + w * y),     0,
			2 * (x * y + w * z),     1 - 2 * (x * x + z * z), 2 * (y * z - w * x),     0,
			2 * (x * z - w * y),     2 * (y * z + w * x),     1 - 2 * (x * x + y * y), 0,
			0,                       0,                       0,                       1
		};
	}

	else { // R == 3 && C == 3
		return {
			1 - 2 * (y * y + z * z), 2 * (x * y - w * z),     2 * (x * z + w * y),
			2 * (x * y + w * z),     1 - 2 * (x * x + z * z), 2 * (y * z - w * x),
			2 * (x * z - w * y),     2 * (y * z + w * x),     1 - 2 * (x * x + y * y)
		};
	}
}

// assumes radians
auto EulerAngles2RotationMatrix(VectorType<3> auto const& crEuler) -> Matrix<ValueType<decltype(crEuler)>, 4, 4> {
	using T = ValueType<decltype(crEuler)>;
	Matrix<T, 4, 4> result = Matrix<T, 4, 4>::Identity();

	LYRotate(crEuler[1], result);
	LXRotate(crEuler[0], result);
	LZRotate(crEuler[2], result);

	return result;
}

auto Homogeneous(VectorType<3> auto const& crVec) -> Vector<ValueType<decltype(crVec)>, 4> {
	using T = ValueType<decltype(crVec)>;
	return Vector<T, 4>(crVec, static_cast<T>(1));
}

auto LookAtMatrix(VectorType<3> auto const& crEye, VectorType<3> auto const& crTarget,
	VectorType<3> auto const& crWorldUp) -> Matrix<ValueType<decltype(crEye)>, 4, 4> {
	using T = ValueType<decltype(crEye)>;

	Vector<T, 3> localz = Normalize(crEye - crTarget);
	Vector<T, 3> localx = Normalize(Cross(crWorldUp, localz));
	Vector<T, 3> localy = Cross(localz, localx);

	Matrix<T, 4, 4> translation = TranslationMatrix(-crEye[0], -crEye[1], -crEye[2]);
	Matrix<T, 4, 4> rotation = {
		localx[0], localx[1], localx[2], 0,
		localy[0], localy[1], localy[2], 0,
		localz[0], localz[1], localz[2], 0,
		0,         0,         0,         1
	};

	return Mul(rotation, translation);
}

template <Axis Dest, Axis Source>
auto ShearMatrix(Numeric auto const cFactor) -> Matrix<std::remove_cvref_t<decltype(cFactor)>, 4, 4> {
	using T = std::remove_cvref_t<decltype(cFactor)>;

	Matrix<T, 4, 4> result = Matrix<T, 4, 4>::Identity();
	result(Dest, Source) = cFactor;

	return result;
}

#endif