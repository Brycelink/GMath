#ifndef MATRIX_GMATH_H
#define MATRIX_GMATH_H

#include <initializer_list>
#include <sstream>
#include <iostream>

#include <cmath>

#include <xmmintrin.h>

#include "GMath/Vector.h"

#include "GMath/GMath_Defs.h"

//USING COLUMN MAJOR FOR MATRICES -- OPENGL COMPATIBLE
namespace GMath{

    static const bool s_ColumnMajor = true;

    template<typename _T, uint16_t rows, uint16_t cols, bool cMajor>
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

        MatrixType(_T value)
        {
            for(uint16_t i = 0; i < rows*cols; i++)
                elements[i] = value;
        }

        MatrixType(_T data[rows*cols]){
            memcpy(elements, data, rows*cols*sizeof(_T));
        }

        template<bool isCMajor>
        MatrixType(const MatrixType<_T, rows, cols, isCMajor>& other)
        {
            if(isCMajor == cMajor)
                memcpy(elements, other.elements, rows*cols*sizeof(_T));
            else{
                if(isCMajor)
                {
                    for(uint16_t i = 0; i < rows; i++){
                        for(uint16_t j = 0; j < cols; j++){
                            elements[j + i * cols] = other.elements[i + j * rows];
                        }
                    }
                }else{
                    for(uint16_t i = 0; i < rows; i++){
                        for(uint16_t j = 0; j < cols; j++){
                            elements[i + j * rows] = other.elements[j + i * cols];
                        }
                    }
                }

            }
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

        static MatrixType<_T, rows, cols, cMajor> Diagonal(_T value)
        {
            MatrixType<_T, rows, cols, cMajor> mat;
            for(uint16_t i = 0; i < rows; i++){
                for(uint16_t j = 0; j < cols; j++){
                    if(i == j)
                        mat.elements[j + i * rows] = value;
                    else
                        mat.elements[j + i * rows] = 0;
                }
            }
            return mat;
        }

        _T At(uint32_t row, uint32_t col)
        {
            if(cMajor)
                return elements[row + col * rows];
            else
                return elements[col + row * cols];
        }

        MatrixType<_T, rows, cols, cMajor>& Add(const MatrixType<_T, rows, cols, cMajor>& other)
        {
            for(int i = 0; i < rows*cols; i++)
                elements[i] += other.elements[i];
            return (*this);
        }

        MatrixType<_T, rows, cols, cMajor>& Subtract(const MatrixType<_T, rows, cols, cMajor>& other)
        {
            for(int i = 0; i < rows*cols; i++)
                elements[i] -= other.elements[i];
            return (*this);
        }

        MatrixType<_T, cols, rows, cMajor> Transpose()
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
        MatrixType<_T, rows, n, cMajor> Multiply(const MatrixType<_T, cols, n, cMajor>& other) const
        {
            _T data[rows*n];
            if(cMajor)
            {
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
            }else{
                for (uint16_t y = 0; y < rows; y++)
                {
                    for(uint16_t x = 0; x < n; x++)
                    {
                        _T sum = 0.0f;
                        for(uint16_t e= 0; e < cols; e++)
                        {
                            sum += elements[e + y * cols] * other.elements[x + e * n];
                        }
                        data[x + y * n] = sum;
                    }
                }
            }

            return MatrixType<_T, rows, n, cMajor>(data);
        }

        friend MatrixType<_T, rows, cols, cMajor> operator+(MatrixType<_T, rows, cols, cMajor> left, const MatrixType<_T, rows, cols, cMajor>& right)
        {
            return left.Add(right);
        }

        friend MatrixType<_T, rows, cols, cMajor> operator-(MatrixType<_T, rows, cols, cMajor> left, const MatrixType<_T, rows, cols, cMajor>& right)
        {
            return left.Subtract(right);
        }

        friend MatrixType<_T, rows, cols, cMajor>& operator+=(MatrixType<_T, rows, cols, cMajor>& left, const MatrixType<_T, rows, cols, cMajor>& right)
        {
            return left.Add(right);
        }

        friend MatrixType<_T, rows, cols, cMajor>& operator-=(MatrixType<_T, rows, cols, cMajor>& left, const MatrixType<_T, rows, cols, cMajor>& right)
        {
            return left.Subtract(right);
        }

        friend bool operator==(const MatrixType<_T, rows, cols, cMajor>& left, const MatrixType<_T, rows, cols, cMajor>& right)
        {
            for(int i = 0; i < rows * cols; i++)
            {
                if(left.elements[i] != right.elements[i])
                    return false;
            }
            return true;
        }

        friend bool operator!=(const MatrixType<_T, rows, cols, cMajor>& left, const MatrixType<_T, rows, cols, cMajor>& right)
        {
            return !(left == right);
        }

        void Print()
        {
            if(cMajor)
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
            }else{
                for(int i = 0; i < rows; i++)
                {
                    printf("[\t");
                    for(int j = 0; j < cols; j++)
                    {
                        printf("%f\t", elements[j + i * cols]);
                    }
                    printf("]\n");
                }
            }
        }

        friend std::ostream& operator<<(std::ostream& stream, const MatrixType<_T, rows, cols, cMajor>& mat)
        {
            if(cMajor)
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
            }else{
                for(int i = 0; i < rows; i++)
                {
                    stream << "[\t";
                    for(int j = 0; j < cols; j++)
                    {
                        stream<< mat.elements[j + i * cols] << "\t";
                    }
                    stream << "]\n";
                }
                return stream;
            }

        }

        static MatrixType<_T, 4, 4, cMajor> Orthographic(const float left, const float right, const float bottom, const float top, const float n, const float f)
        {
            MatrixType<_T, 4, 4, cMajor> result = GMath::MatrixType<_T, 4, 4, cMajor>::Diagonal(1.0f);

            result.elements[0 + 0 * 4] = (_T)(2.0f / (right - left));
            result.elements[1 + 1 * 4] = (_T)(2.0f / (top - bottom));
            result.elements[2 + 2 * 4] = (_T)(-2.0f / (f - n));

            if(cMajor){
                result.elements[0 + 3 * 4] = (_T)(-(right + left)/(right - left));
                result.elements[1 + 3 * 4] = (_T)(-(top + bottom)/(top - bottom));
                result.elements[2 + 3 * 4] = (_T)(-(f + n)/(f - n));
            }else{
                result.elements[3 + 0 * 4] = (_T)(-(right + left)/(right - left));
                result.elements[3 + 1 * 4] = (_T)(-(top + bottom)/(top - bottom));
                result.elements[3 + 2 * 4] = (_T)(-(f + n)/(f - n));
            }

            return result;
        }

        static MatrixType<_T, 4, 4, cMajor> Perspective(const float fov, const float aspectRatio, const float n, const float f)
        {
            MatrixType<_T, 4, 4, cMajor> result = GMath::MatrixType<_T, 4, 4, cMajor>::Diagonal(1.0f);

            float q = 1.0f / (float)tan(toRadians(fov / 2.0f));
            float a = q / aspectRatio;

            float b = (n + f) / (n - f);
            float c = (2.0f * n * f) / (n - f);

            result.elements[0 + 0 * 4] = a;
            result.elements[1 + 1 * 4] = q;
            result.elements[2 + 2 * 4] = b;

            if(cMajor){
                result.elements[3 + 2 *4] = -1.0f;
                result.elements[2 + 3 *4] = c;
            }else{
                result.elements[2 + 3 *4] = -1.0f;
                result.elements[3 + 2 *4] = c;
            }

            return result;
        }

        static MatrixType<_T, 4, 4, cMajor> Translation(const VectorType<_T, 3> &translation)
        {
            MatrixType<_T, 4, 4, cMajor> result = GMath::MatrixType<_T, 4, 4, cMajor>::Diagonal(1.0f);
            if(cMajor){
                result.elements[12] = translation.elements[0];
                result.elements[13] = translation.elements[1];
                result.elements[14] = translation.elements[2];
            }else{
                result.elements[3] = translation.elements[0];
                result.elements[7] = translation.elements[1];
                result.elements[11] = translation.elements[2];
            }

            return result;
        }

        static MatrixType<_T, 4, 4, cMajor> Scale(const VectorType<_T, 3> &scale)
        {
            MatrixType<_T, 4, 4, cMajor> result = GMath::MatrixType<_T, 4, 4, cMajor>::Diagonal(1.0f);

            result.elements[0] = scale.elements[0];
            result.elements[5] = scale.elements[1];
            result.elements[10] = scale.elements[2];

            return result;
        }

        static MatrixType<_T, 4, 4, cMajor> Rotation(const float& angle, const VectorType<_T, 3> &axis)
        {
            MatrixType<_T, 4, 4, cMajor> result = GMath::MatrixType<_T, 4, 4, cMajor>::Diagonal(1.0f);

            float r = toRadians(angle);
            float c = cos(r);
            float s = sin(r);
            float omc = 1.0f - c;

            float x = axis.elements[0];
            float y = axis.elements[1];
            float z = axis.elements[2];

            if(cMajor)
            {
                result.elements[0] = x * omc + c;
                result.elements[1] = y * x * omc + z * s;
                result.elements[2] = x * z * omc - y * s;

                result.elements[4] = x * y * omc - z * s;
                result.elements[5] = y * omc + c;
                result.elements[6] = y * z * omc + x * s;

                result.elements[8] = x * z * omc + y * s;
                result.elements[9] = y * z * omc - x * s;
                result.elements[10] = z * omc + c;
            }else{
                result.elements[0] = x * omc + c;
                result.elements[4] = y * x * omc + z * s;
                result.elements[8] = x * z * omc - y * s;

                result.elements[1] = x * y * omc - z * s;
                result.elements[5] = y * omc + c;
                result.elements[9] = y * z * omc + x * s;

                result.elements[2] = x * z * omc + y * s;
                result.elements[6] = y * z * omc - x * s;
                result.elements[10] = z * omc + c;
            }

            return result;
        }
    };

    template<typename _T, uint16_t rows, uint16_t cols, bool cMajor>
    static MatrixType<_T, cols, rows, cMajor> Transpose(const MatrixType<_T, cols, rows, cMajor>& matrix)
    {
        _T data[rows*cols];
        for(uint16_t i = 0; i < rows; i++){
            for(uint16_t j = 0; j < cols; j++){
                data[j + i * cols] = matrix.elements[i + j * rows];
            }
        }
        return MatrixType<_T, cols, rows, cMajor>(data);
    }

    template<typename _T, bool cMajor>
    static MatrixType<_T, 4, 4, cMajor> Multiply(const MatrixType<_T, 4, 4, cMajor>& matrix, const MatrixType<_T, 4, 4, cMajor>& other)
    {
        _T res[16];
        if(cMajor){
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
        }else{
            __m128 rowA1 = _mm_loadu_ps(&other.elements[0]);
            __m128 rowA2 = _mm_loadu_ps(&other.elements[4]);
            __m128 rowA3 = _mm_loadu_ps(&other.elements[8]);
            __m128 rowA4 = _mm_loadu_ps(&other.elements[12]);

            for(int i=0; i<4; i++) {
                __m128 brod1 = _mm_set1_ps(matrix.elements[4*i + 0]);
                __m128 brod2 = _mm_set1_ps(matrix.elements[4*i + 1]);
                __m128 brod3 = _mm_set1_ps(matrix.elements[4*i + 2]);
                __m128 brod4 = _mm_set1_ps(matrix.elements[4*i + 3]);
                __m128 row = _mm_add_ps(
                            _mm_add_ps(
                                _mm_mul_ps(brod1, rowA1),
                                _mm_mul_ps(brod2, rowA2)),
                            _mm_add_ps(
                                _mm_mul_ps(brod3, rowA3),
                                _mm_mul_ps(brod4, rowA4)));
                _mm_store_ps(&res[4*i], row);
            }
        }

        return res;
    }

    template<typename _T, bool cMajor>
    static MatrixType<_T, 4, 4, cMajor> operator*(const MatrixType<_T, 4, 4, cMajor>& matrix, const MatrixType<_T, 4, 4, cMajor>& other)
    {
        return Multiply(matrix, other);
    }

    template<typename _T, uint16_t rows, uint16_t cols, uint16_t n, bool cMajor>
    static MatrixType<_T, rows, n, cMajor> operator*(const MatrixType<_T, rows, cols, cMajor>& left, const MatrixType< _T, cols, n, cMajor>& right)
    {
        return left.Multiply(right);
    }

template<typename _T, uint16_t rows, uint16_t cols, bool cMajor>
    VectorType<_T, rows> Multiply(const MatrixType<_T, rows, cols, cMajor>& matrix, const VectorType<_T, cols>& vec)
    {
        _T data[rows];
        if(cMajor){
            for (int y = 0; y < rows; y++)
            {
                _T sum = 0.0f;
                for(int e= 0; e < cols; e++)
                {
                    sum += matrix.elements[y + e * rows] * vec.elements[e];
                }
                data[y] = sum;
            }
        }else{
            for (int y = 0; y < rows; y++)
            {
                _T sum = 0.0f;
                for(int e= 0; e < cols; e++)
                {
                    sum += matrix.elements[e + y * cols] * vec.elements[e];
                }
                data[y] = sum;
            }
        }

        return VectorType<_T, rows>(data);
    }

    template<typename _T, uint16_t rows, uint16_t cols, bool cMajor>
    VectorType<_T, rows> operator*(const MatrixType<_T, rows, cols, cMajor>& matrix, const VectorType<_T, cols>& vec){
        return Multiply(matrix, vec);
    }

    template<typename _T, uint16_t n, bool cMajor>
    static MatrixType<_T, n-1, n-1, cMajor> Minor(const MatrixType<_T, n, n, cMajor>& matrix, uint16_t row, uint16_t col)
    {
        MatrixType<_T, n-1, n-1, cMajor> result;
        uint16_t pos = 0;

        if(cMajor){
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
        }else{
            for(uint16_t j = 0; j < n; j++)
            {
                for (uint16_t i  = 0; i < n; i++)
                {
                    if ((row != j) && (col != i))
                    {
                        result.elements[pos] = matrix.elements[i + j * n];
                        pos++;
                    }
                }
            }
        }

        return result;
    }

    template<typename _T, uint16_t n, bool cMajor>
    static GMATH_TYPE_PRECISION Determinant(const MatrixType<_T, n, n, cMajor>& matrix)
    {
        GMATH_TYPE_PRECISION det = 0;
        uint16_t k = 1;

        if(cMajor){
            for(uint16_t j = 0; j < n; j++)
                det += matrix.elements[k + j * n ] * Cofactor(matrix, k, j);
        }else{
            for(uint16_t j = 0; j < n; j++)
                det += matrix.elements[j + k * n ] * Cofactor(matrix, k, j);
        }
        return det;
    }

    template<bool cMajor> inline GMATH_TYPE_PRECISION Determinant(const MatrixType<GMATH_TYPE_PRECISION, 2, 2, cMajor>& matrix)
    {
        return (matrix.elements[0]*matrix.elements[3]) - (matrix.elements[1]*matrix.elements[2]);
    }

    template<bool cMajor> inline GMATH_TYPE_PRECISION Determinant(const MatrixType<GMATH_TYPE_PRECISION, 1, 1, cMajor>& matrix)
    {
        return matrix.elements[0];
    }

    template<typename _T, uint16_t n, bool cMajor>
    static GMATH_TYPE_PRECISION Cofactor(const MatrixType<_T, n, n, cMajor>& matrix, uint16_t row, uint16_t col)
    {
        GMATH_TYPE_PRECISION value;
        int convention;

        if((row+col)%2 == 0)
            convention = 1;
        else
            convention = -1;

        MatrixType<_T, n-1, n-1, cMajor> midPoint = Minor(matrix, row, col);
        value = convention * Determinant(midPoint);

        return value;
    }

    template<typename _T, uint16_t n, bool cMajor>
    static bool IsInvertible(MatrixType<_T, n, n, cMajor>& matrix){
        if(Determinant(matrix) == 0)
            return false;
        return true;
    }

    template<typename _T, uint16_t n, bool cMajor>
    static MatrixType<_T, n, n, cMajor> InverseAdj(const MatrixType<_T, n, n, cMajor>& matrix)
    {
        GMATH_TYPE_PRECISION determinant = Determinant(matrix);
        if(determinant == 0){
            printf("Matrix is Invertible!\n");
            return matrix;
        }
        GMATH_TYPE_PRECISION invDet = 1/determinant;
        _T data[n*n];
        if(cMajor){
            for(uint16_t i = 0; i < n; i++){
                for(uint16_t j = 0; j < n; j++){
                    data[i + j * n] = Cofactor(matrix, j, i) * invDet;
                }
            }
        }else{
            for(uint16_t i = 0; i < n; i++){
                for(uint16_t j = 0; j < n; j++){
                    data[j + i * n] = Cofactor(matrix, j, i) * invDet;
                }
            }
        }

        return MatrixType<_T, n, n, cMajor>(data);
    }

    template<typename _T, uint16_t n, bool cMajor>
    static MatrixType<_T, n, n, cMajor> Inverse(MatrixType<_T, n, n, cMajor>& matrix)
    {
        if(Determinant(matrix) == 0){
            printf("Matrix is Invertible!\n");
            std::cout<<"Here"<<std::endl;
            return matrix;
        }

        _T data[n*n*2];
        MatrixType<_T, n, n, cMajor> I = GMath::MatrixType<_T, n, n, cMajor>::Diagonal(1);

        for(uint16_t i = 0; i < n*n*2; i++)
        {
            if(i<n*n){
                data[i] = matrix.elements[i];
            }else
                data[i] = I.elements[i-n*n];
        }

        for (uint16_t j = 0; j < n; j++){
            _T largest = fabs(data[j + j * n]);
            uint16_t index = j;

            for (uint16_t i = j; i < n; i++){
                if ((fabs(data[i + j * n]) > largest) && (fabs(data[i + j * n]) != 0)){
                    largest = fabs(data[i + j * n]);
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
                        for(uint16_t r = j; r < n*2; r++){
                            data[i + r * n] += refer * data[j + r * n];
                        }
                    }
                }
            }
        }
        _T res[n*n];
        for(uint16_t i = n*n; i < n*n*2; i++)
            res[i-n*n] = data[i];
        return MatrixType<_T, n, n, cMajor>(res);
    }

    template<typename _T, uint16_t rows, uint16_t cols>
    using MatrixTypeC = MatrixType<GMATH_TYPE_PRECISION, rows, cols, true>;

    template<typename _T, uint16_t rows, uint16_t cols>
    using MatrixTypeR = MatrixType<GMATH_TYPE_PRECISION, rows, cols, false>;

    template<typename _T, uint16_t rows, uint16_t cols, bool cMajor>
    const _T* value_ptr(const MatrixType< _T, rows, cols, cMajor>& matrix) {return &matrix.elements[0];}

}

#endif // MATRIX_GMATH_H

