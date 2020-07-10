#ifndef GMATH_GMATH_H
#define GMATH_GMATH_H

#include "GMath/Matrix.h"
#include "GMath/Vector.h"
#include "GMath/Quaternion.h"

namespace GMath
{

    template<typename _T, uint16_t rows, uint16_t cols>
    VectorType<_T, rows> Multiply(const MatrixType<_T, rows, cols>& matrix, const VectorType<_T, cols>& vec)
    {
        _T data[rows];
        for (int y = 0; y < rows; y++)
        {
            for(int x = 0; x < 1; x++)
            {
                _T sum = 0.0f;
                for(int e= 0; e < cols; e++)
                {
                    sum += matrix.elements[y + e * rows] * vec.elements[e + x * cols];
                }
                data[y + x * rows] = sum;
            }
        }
        return VectorType<_T, rows>(data);
    }

    template<typename _T, uint16_t rows, uint16_t cols>
    VectorType<_T, rows> operator*(const MatrixType<_T, rows, cols>& matrix, const VectorType<_T, cols>& vec){
        return Multiply(matrix, vec);
    }

    template<typename _T = GMATH_TYPE_PRECISION>
    static MatrixType<_T, 4, 4> Translation(const VectorType<_T, 3> &translation)
    {
        MatrixType<_T, 4, 4> result(1.0f);

        result.elements[0 + 3 * 4] = translation.elements[0];
        result.elements[1 + 3 * 4] = translation.elements[1];
        result.elements[2 + 3 * 4] = translation.elements[2];

        return result;
    }

    template<typename _T = GMATH_TYPE_PRECISION>
    static MatrixType<_T, 4, 4> Scale(const VectorType<_T, 3> &scale)
    {
        MatrixType<_T, 4, 4> result(1.0f);

        result.elements[0 + 0 * 4] = scale.elements[0];
        result.elements[1 + 1 * 4] = scale.elements[1];
        result.elements[2 + 2 * 4] = scale.elements[2];

        return result;
    }

    template<typename _T = GMATH_TYPE_PRECISION>
    static MatrixType<_T, 4, 4> Rotation(const float& angle, const VectorType<_T, 3> &axis)
    {
        MatrixType<_T, 4, 4> result(1.0f);

        float r = toRadians(angle);
        float c = cos(r);
        float s = sin(r);
        float omc = 1.0f - c;

        float x = axis.elements[0];
        float y = axis.elements[1];
        float z = axis.elements[2];

        result.elements[0 + 0 * 4] = x * omc + c;
        result.elements[1 + 0 * 4] = y * x * omc + z * s;
        result.elements[2 + 0 * 4] = x * z * omc - y * s;

        result.elements[0 + 1 * 4] = x * y * omc - z * s;
        result.elements[1 + 1 * 4] = y * omc + c;
        result.elements[2 + 1 * 4] = y * z * omc + x * s;

        result.elements[0 + 2 * 4] = x * z * omc + y * s;
        result.elements[1 + 2 * 4] = y * z * omc - x * s;
        result.elements[2 + 2 * 4] = z * omc + c;

        return result;
    }
}

#endif // GMATH_GMATH_H
