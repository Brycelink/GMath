#ifndef MATRIX_GMATH_H
#define MATRIX_GMATH_H

#include <initializer_list>
#include <sstream>
#include <iostream>

#include <cmath>

#include <xmmintrin.h>

#include "Logs.h"

//USING COLUMN MAJOR FOR MATRICES -- OPENGL COMPATIBLE
namespace GMath{

    template<typename _T, uint16_t rows, uint16_t cols>
    struct MatrixType
    {
    public:
        _T elements[rows*cols];

    public:
        MatrixType()
        {
            for(uint16_t i = 0; i < rows*cols; i++)
                elements[i] = 0.0f;
        }

        MatrixType(_T diag, bool fill = false)
        {
            if(fill){
                for(uint16_t i = 0; i < rows*cols; i++)
                    elements[i] = diag;
                return;
            }
            for(uint16_t i = 0; i < rows; i++){
                for(uint16_t j = 0; j < cols; j++){
                    if(i == j)
                        elements[j + i * rows] = diag;
                    else
                        elements[j + i * rows] = 0;
                }
            }
        }

        MatrixType(_T data[rows*cols]){
            memcpy(elements, data, rows*cols*sizeof(_T));
        }

        MatrixType(const MatrixType<_T, rows, cols>& other)
        {
            memcpy(elements, other.elements, rows*cols*sizeof(_T));
        }

        MatrixType(const std::initializer_list<_T>& data){
            _T* temp = (_T*)data.begin();
            for(uint16_t i = 0; i < data.size(); i++){
                elements[i] = temp[i];
            }
        }

        MatrixType(const std::initializer_list<std::initializer_list<_T>>& data)
        {
            int i = 0, j = 0;
            for (const auto& l : data)
            {
                for (const auto& v : l)
                {
                   elements[j + i * rows] = v;
                   ++j;
                }
                j = 0;
                ++i;
            }
        }

        _T At(uint32_t row, uint32_t col)
        {
            return elements[row + col * rows];
        }

        MatrixType<_T, rows, cols>& Add(const MatrixType<_T, rows, cols>& other)
        {
            for(int i = 0; i < rows*cols; i++)
                elements[i] += other.elements[i];
            return (*this);
        }

        MatrixType<_T, rows, cols>& Subtract(const MatrixType<_T, rows, cols>& other)
        {
            for(int i = 0; i < rows*cols; i++)
                elements[i] -= other.elements[i];
            return (*this);
        }

        MatrixType<_T, cols, rows> Transpose()
        {
            _T data[rows*cols];
            for(uint16_t i = 0; i < rows; i++){
                for(uint16_t j = 0; j < cols; j++){
                    data[j + i * cols] = elements[i + j * rows];
                }
            }
            (*this) = data;
            return (*this);
        }

        template<uint16_t n>
        MatrixType<_T, rows, n> Multiply(const MatrixType<_T, cols, n>& other) const
        {
            _T data[rows*n];
            for (uint16_t y = 0; y < rows; y++)
            {
                for(uint16_t x = 0; x < n; x++)
                {
                    _T sum = 0.0f;
                    for(uint16_t e= 0; e < cols; e++)
                    {
                        sum += elements[y + e * rows] * other.elements[e + x * cols];
                    }
                    data[y + x * rows] = sum;
                }
            }
            return MatrixType<_T, rows, n>(data);
        }

        friend MatrixType<_T, rows, cols> operator+(MatrixType<_T, rows, cols> left, const MatrixType<_T, rows, cols>& right)
        {
            return left.Add(right);
        }

        friend MatrixType<_T, rows, cols> operator-(MatrixType<_T, rows, cols> left, const MatrixType<_T, rows, cols>& right)
        {
            return left.Subtract(right);
        }

        friend MatrixType<_T, rows, cols>& operator+=(MatrixType<_T, rows, cols>& left, const MatrixType<_T, rows, cols>& right)
        {
            return left.Add(right);
        }

        friend MatrixType<_T, rows, cols>& operator-=(MatrixType<_T, rows, cols>& left, const MatrixType<_T, rows, cols>& right)
        {
            return left.Subtract(right);
        }

        friend bool operator==(const MatrixType<_T, rows, cols>& left, const MatrixType<_T, rows, cols>& right)
        {
            for(int i = 0; i < rows * cols; i++)
            {
                if(left.elements[i] != right.elements[i])
                    return false;
            }
            return true;
        }

        friend bool operator!=(const MatrixType<_T, rows, cols>& left, const MatrixType<_T, rows, cols>& right)
        {
            return !(left == right);
        }

        void Print()
        {
            for(int i = 0; i < rows; i++)
            {
                printf("[\t");
                for(int j = 0; j < cols; j++)
                {
                    printf("%f\t", elements[i + j * rows]);
                }
                printf("]\n");
            }
        }

        friend std::ostream& operator<<(std::ostream& stream, const MatrixType<_T, rows, cols>& mat)
        {
            for(int i = 0; i < rows; i++)
            {
                stream << "[\t";
                for(int j = 0; j < cols; j++)
                {
                    stream<< mat.elements[i + j * rows] << "\t";
                }
                stream << "]\n";
            }
            return stream;
        }
    };

    template<typename _T, uint16_t rows, uint16_t cols>
    static MatrixType<_T, cols, rows> Transpose(const MatrixType<_T, cols, rows>& matrix)
    {
        _T data[rows*cols];
        for(uint16_t i = 0; i < rows; i++){
            for(uint16_t j = 0; j < cols; j++){
                data[j + i * cols] = matrix.elements[i + j * rows];
            }
        }
        return MatrixType<_T, cols, rows>(data);
    }

    template<typename _T>
    static MatrixType<_T, 4, 4> Multiply(const MatrixType<_T, 4, 4>& matrix, const MatrixType<_T, 4, 4>& other)
    {
        _T res[16];
        __m128 rowA1 = _mm_loadu_ps(&matrix.elements[0]);
        __m128 rowA2 = _mm_loadu_ps(&matrix.elements[4]);
        __m128 rowA3 = _mm_loadu_ps(&matrix.elements[8]);
        __m128 rowA4 = _mm_loadu_ps(&matrix.elements[12]);

        for(int i=0; i<4; i++) {
            __m128 brod1 = _mm_set1_ps(other.elements[4*i + 0]);
            __m128 brod2 = _mm_set1_ps(other.elements[4*i + 1]);
            __m128 brod3 = _mm_set1_ps(other.elements[4*i + 2]);
            __m128 brod4 = _mm_set1_ps(other.elements[4*i + 3]);
            __m128 row = _mm_add_ps(
                        _mm_add_ps(
                            _mm_mul_ps(brod1, rowA1),
                            _mm_mul_ps(brod2, rowA2)),
                        _mm_add_ps(
                            _mm_mul_ps(brod3, rowA3),
                            _mm_mul_ps(brod4, rowA4)));
            _mm_store_ps(&res[4*i], row);
        }
        return res;
    }

    template<typename _T>
    static MatrixType<_T, 4, 4> operator*(const MatrixType<_T, 4, 4>& matrix, const MatrixType<_T, 4, 4>& other)
    {
        return Multiply(matrix, other);
    }

    template<typename _T, uint16_t rows, uint16_t cols, uint16_t n>
    static MatrixType<_T, rows, n> operator*(const MatrixType<_T, rows, cols>& left, const MatrixType< _T, cols, n>& right)
    {
        return left.Multiply(right);
    }

    template<typename _T, uint16_t n>
    static MatrixType<_T, n-1, n-1> Minor(MatrixType<_T, n, n>& matrix, uint16_t row, uint16_t col)
    {
        MatrixType<_T, n-1, n-1> result;
        uint16_t pos = 0;

        for(uint16_t j = 0; j < n; j++)
        {
            for (uint16_t i  = 0; i < n; i++)
            {
                if ((row != i) && (col != j))
                {
                    result.elements[pos] = matrix.elements[i + j * n];
                    pos++;
                }
            }
        }
        return result;
    }

    template<typename _T, uint16_t n>
    static GMATH_TYPE_PRECISION Determinant(MatrixType<_T, n, n>& matrix)
    {
        GMATH_TYPE_PRECISION det;
        uint16_t k = 1;

        for(uint16_t j = 0; j < n; j++)
            det += matrix.elements[k + j * n] * Cofactor(matrix, k, j);

        return det;
    }

    template<typename _T, uint16_t n>
    static GMATH_TYPE_PRECISION Cofactor(MatrixType<_T, n, n>& matrix, uint16_t row, uint16_t col)
    {
        GMATH_TYPE_PRECISION value;
        int convention;

        if((row+col)%2 == 0)
            convention = 1;
        else
            convention = -1;

        MatrixType<_T, n-1, n-1> midPoint = Minor(matrix, row, col);
        value = convention * Determinant(midPoint);

        return value;
    }

    template<> inline GMATH_TYPE_PRECISION Determinant(MatrixType<GMATH_TYPE_PRECISION, 2, 2>& matrix)
    {
        return (matrix.elements[0]*matrix.elements[3]) - (matrix.elements[1]*matrix.elements[2]);
    }

    template<> inline GMATH_TYPE_PRECISION Determinant(MatrixType<GMATH_TYPE_PRECISION, 1, 1>& matrix)
    {
        return matrix.elements[0];
    }

    template<typename _T, uint16_t n>
    static bool IsInvertible(MatrixType<_T, n, n>& matrix){
        if(Determinant(matrix) == 0)
            return true;
        return false;
    }

    template<typename _T, uint16_t n>
    static MatrixType<_T, n, n> InverseAdj(MatrixType<_T, n, n>& matrix)
    {
        GMATH_TYPE_PRECISION determinant = Determinant(matrix);
        if(determinant == 0){
            printf("Matrix is Invertible!\n");
            return matrix;
        }
        GMATH_TYPE_PRECISION det = 1/determinant;
        _T data[n*n];
        for(uint16_t i = 0; i < n; i++){
            for(uint16_t j = 0; j < n; j++){
                data[i + j * n] = Cofactor(matrix, j, i) * det;
            }
        }
        return MatrixType<_T, n, n>(data);
    }

    template<typename _T, uint16_t n>
    static MatrixType<_T, n, n> Inverse(MatrixType<_T, n, n>& matrix)
    {
        if(Determinant(matrix) == 0){
            printf("Matrix is Invertible!\n");
            return matrix;
        }

        _T data[n*n*2];
        MatrixType<_T, n, n> I(1, true);
        for(uint16_t i = 0; i < n*n*2; i++)
        {
            if(i<n*n){
                data[i] = matrix.elements[i];
            }else
                data[i] = I.elements[i-n*n];
        }

        for (uint16_t j = 0; j < n; j++){
            _T largest = abs(data[j + j * n]);
            uint16_t index = j;

            for (uint16_t i = j; i < n; i++){
                if ((abs(data[i + j * n]) > largest) && (abs(data[i + j * n]) != 0)){
                    largest = abs(data[i + j * n]);
                    index = i;
                }
            }
            if (largest == 0){
                printf("Matrix is Invertible!\n");
                return matrix;
            }else{
                if( index != j){
                    for (uint16_t i = j; i < n*2; i++){
                        _T temp = data[j + i * n];
                        data[j + i * n] = data[index + i * n];
                        data[index + i * n] = temp;
                    }
                }

                _T mainValue = (1/(data[j + j * n]));

                for (uint16_t i = 0; i < n*2; i++)
                    data[j + i * n] *= mainValue;

                for (uint16_t i = 0; i < n; i++){
                    if(i != j){
                        _T refer = -1 * data[i + j * n];
                        for(uint16_t r = 0; r < n*2; r++){
                            data[i + r * n] += refer * data[j + r * n];
                        }
                    }
                }
            }
        }
        _T res[n*n];
        for(uint16_t i = n*n; i < n*n*2; i++)
            res[i-n*n] = data[i];
        return MatrixType<_T, n, n>(res);
    }

    template<typename _T = GMATH_TYPE_PRECISION>
    static MatrixType<_T, 4, 4> Orthographic(const float left, const float right, const float bottom, const float top, const float n, const float f)
    {
        MatrixType<_T, 4, 4> result(1.0f);

        result.elements[0 + 0 * 4] = (_T)(2.0f / (right - left));
        result.elements[1 + 1 * 4] = (_T)(2.0f / (top - bottom));
        result.elements[2 + 2 * 4] = (_T)(-2.0f / (f - n));

        result.elements[0 + 3 * 4] = (_T)(-(right + left)/(right - left));
        result.elements[1 + 3 * 4] = (_T)(-(top + bottom)/(top - bottom));
        result.elements[2 + 3 * 4] = (_T)(-(f + n)/(f - n));

        return result;
    }

    template<typename _T = GMATH_TYPE_PRECISION>
    static MatrixType<_T, 4, 4> Perspective(const float fov, const float aspectRatio, const float n, const float f)
    {
        MatrixType<_T, 4, 4> result(1.0f);

        float q = 1.0f / (float)tan(toRadians(fov / 2.0f));
        float a = q / aspectRatio;

        float b = (n + f) / (n - f);
        float c = (2.0f * n * f) / (n - f);

        result.elements[0 + 0 * 4] = a;
        result.elements[1 + 1 * 4] = q;
        result.elements[2 + 2 * 4] = b;
        result.elements[3 + 2 *4] = -1.0f;
        result.elements[2 + 3 *4] = c;

        return result;
    }


    template<uint16_t rows, uint16_t cols>
    using Matrix = MatrixType<GMATH_TYPE_PRECISION, rows, cols>;

    template<uint16_t n>
    using MatrixSq = MatrixType<GMATH_TYPE_PRECISION, n, n>;

    typedef MatrixSq<4> Matrix4;
    typedef MatrixSq<3> Matrix3;
    typedef MatrixSq<2> Matrix2;
    typedef MatrixSq<1> Matrix1;

    template<typename _T, uint16_t rows, uint16_t cols>
    _T* value_ptr(MatrixType< _T, rows, cols> matrix) {return &matrix.elements[0];}

}

#endif // MATRIX_GMATH_H

