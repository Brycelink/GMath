#ifndef VECTOR_GMATH_H
#define VECTOR_GMATH_H

#include "Logs.h"
#include <cmath>

namespace GMath{

    template<typename _T, uint16_t n>
    struct VectorType
    {
    public:
        VectorType()
        {
            for(uint16_t i = 0; i < n; i++)
                elements[i] = 0.0f;
        }

        VectorType(float value)
        {
            for(uint16_t i = 0; i < n; i++)
                elements[i] = value;
        }

        VectorType(float data[n]){
            memcpy(elements, data, n*sizeof(float));
        }

        VectorType(std::initializer_list<_T> data)
        {
            _T* temp = (_T*)data.begin();
            for(uint16_t i = 0; i < data.size(); i++){
                elements[i] = temp[i];
            }
        }

        GMATH_TYPE_PRECISION Magnitude()
        {
            GMATH_TYPE_PRECISION total = 0.0f;
            for(uint16_t i = 0; i < n; i++)
                total += elements[i] * elements[i];
            return sqrt(total);
        }

        GMATH_TYPE_PRECISION Length(){return Magnitude();}

        GMATH_TYPE_PRECISION SqMagnitude()
        {
            GMATH_TYPE_PRECISION total = 0.0f;
            for(uint16_t i = 0; i < n; i++)
                total += elements[i] * elements[i];
            return total;
        }

        GMATH_TYPE_PRECISION SqLength(){return SqMagnitude();}

        VectorType<_T, n>& Normalize()
        {
            GMATH_TYPE_PRECISION length = Magnitude();
            for(uint16_t i = 0; i < n; i++)
                elements[i] /= length;
            return (*this);
        }

        VectorType<_T, n> Negate()
        {
            for(uint16_t i = 0; i < n; i++)
                elements[i] = -elements[i];
            return (*this);
        }

        void Print()
        {
            printf("( ");
            for(uint16_t i = 0; i < n; i++){
                printf("%f", elements[i]);
                if(i != n-1)
                    printf(", ");
            }
            printf(")\n");
        }

    public:
        _T elements[n];
    };

//VECTOR 4
    template <typename _T>
    struct VectorType<_T, 4>
    {
    public:
        VectorType()
        {
            for(uint16_t i = 0; i < 4; i++)
                elements[i] = 0.0f;
        }

        VectorType(float value)
        {
            for(uint16_t i = 0; i < 4; i++)
                elements[i] = value;
        }

        VectorType(float _x, float _y, float _z, float _w)
        {
            x = _x; y = _y;
            z = _z; w = _w;
        }

        VectorType(float data[4]){
            memcpy(elements, data, 4*sizeof(float));
        }

        VectorType(std::initializer_list<_T> data)
        {
            _T* temp = (_T*)data.begin();
            for(uint16_t i = 0; i < data.size(); i++){
                elements[i] = temp[i];
            }
        }

        GMATH_TYPE_PRECISION Magnitude()
        {
            GMATH_TYPE_PRECISION total = 0.0f;
            for(uint16_t i = 0; i < 4; i++)
                total += elements[i] * elements[i];
            return sqrt(total);
        }

        GMATH_TYPE_PRECISION Length(){return Magnitude();}

        GMATH_TYPE_PRECISION SqMagnitude()
        {
            GMATH_TYPE_PRECISION total = 0.0f;
            for(uint16_t i = 0; i < 4; i++)
                total += elements[i] * elements[i];
            return total;
        }

        GMATH_TYPE_PRECISION SqLength(){return SqMagnitude();}

        VectorType<_T, 4>& Normalize()
        {
            GMATH_TYPE_PRECISION length = Magnitude();
            for(uint16_t i = 0; i < 4; i++)
                elements[i] /= length;
            return (*this);
        }

        VectorType<_T, 4> Negate()
        {
            for(uint16_t i = 0; i < 4; i++)
                elements[i] = -elements[i];
            return (*this);
        }

        void Print()
        {
            printf("( ");
            for(uint16_t i = 0; i < 4; i++){
                printf("%f", elements[i]);
                if(i != 3)
                    printf(", ");
            }
            printf(")\n");
        }

    public:
        union{
            struct{ _T x, y, z, w;};
            struct{ _T r, g, b, a;};
            _T elements[4];
        };
    };


//VECTOR 3
    template <typename _T>
    struct VectorType<_T, 3>
    {
    public:
        VectorType()
        {
            for(uint16_t i = 0; i < 3; i++)
                elements[i] = 0.0f;
        }

        VectorType(float value)
        {
            for(uint16_t i = 0; i < 3; i++)
                elements[i] = value;
        }

        VectorType(float data[3]){
            memcpy(elements, data, 3*sizeof(float));
        }

        VectorType(float _x, float _y, float _z)
        {
            x = _x;
            y = _y;
            z = _z;
        }

        VectorType(std::initializer_list<_T> data)
        {
            _T* temp = (_T*)data.begin();
            for(uint16_t i = 0; i < data.size(); i++){
                elements[i] = temp[i];
            }
        }

        GMATH_TYPE_PRECISION Magnitude()
        {
            GMATH_TYPE_PRECISION total = 0.0f;
            for(uint16_t i = 0; i < 3; i++)
                total += elements[i] * elements[i];
            return sqrt(total);
        }

        GMATH_TYPE_PRECISION Length(){return Magnitude();}

        GMATH_TYPE_PRECISION SqMagnitude()
        {
            GMATH_TYPE_PRECISION total = 0.0f;
            for(uint16_t i = 0; i < 3; i++)
                total += elements[i] * elements[i];
            return total;
        }

        GMATH_TYPE_PRECISION SqLength(){return SqMagnitude();}

        VectorType<_T, 3>& Normalize()
        {
            GMATH_TYPE_PRECISION length = Magnitude();
            for(uint16_t i = 0; i < 3; i++)
                elements[i] /= length;
            return (*this);
        }

        VectorType<_T, 3> Negate()
        {
            for(uint16_t i = 0; i < 3; i++)
                elements[i] = -elements[i];
            return (*this);
        }

        void Print()
        {
            printf("( ");
            for(uint16_t i = 0; i < 3; i++){
                printf("%f", elements[i]);
                if(i != 2)
                    printf(", ");
            }
            printf(")\n");
        }
    public:
        union{
            struct{ _T x, y, z;};
            _T elements[3];
        };
    };

//VECTOR 2
    template <typename _T>
    struct VectorType<_T, 2>
    {
    public:
        VectorType()
        {
            for(uint16_t i = 0; i < 2; i++)
                elements[i] = 0.0f;
        }

        VectorType(float value)
        {
            for(uint16_t i = 0; i < 2; i++)
                elements[i] = value;
        }

        VectorType(float data[2]){
            memcpy(elements, data, 2*sizeof(float));
        }
        VectorType(float _x, float _y)
        {
            x = _x;
            y = _y;
        }

        VectorType(std::initializer_list<_T> data)
        {
            _T* temp = (_T*)data.begin();
            for(uint16_t i = 0; i < data.size(); i++){
                elements[i] = temp[i];
            }
        }

        GMATH_TYPE_PRECISION Magnitude()
        {
            GMATH_TYPE_PRECISION total = 0.0f;
            for(uint16_t i = 0; i < 2; i++)
                total += elements[i] * elements[i];
            return sqrt(total);
        }

        GMATH_TYPE_PRECISION Length(){return Magnitude();}

        GMATH_TYPE_PRECISION SqMagnitude()
        {
            GMATH_TYPE_PRECISION total = 0.0f;
            for(uint16_t i = 0; i < 2; i++)
                total += elements[i] * elements[i];
            return total;
        }

        GMATH_TYPE_PRECISION SqLength(){return SqMagnitude();}

        VectorType<_T, 2>& Normalize()
        {
            GMATH_TYPE_PRECISION length = Magnitude();
            for(uint16_t i = 0; i < 2; i++)
                elements[i] /= length;
            return (*this);
        }

        VectorType<_T, 2> Negate()
        {
            for(uint16_t i = 0; i < 2; i++)
                elements[i] = -elements[i];
            return (*this);
        }

        void Print()
        {
            printf("( ");
            for(uint16_t i = 0; i < 2; i++){
                printf("%f", elements[i]);
                if(i != 1)
                    printf(", ");
            }
            printf(")\n");
        }
    public:
        union{
            struct{ _T x, y;};
            _T elements[2];
        };
    };

    template<typename _T, uint16_t n>
    static VectorType<_T, n> Add(const VectorType<_T, n>& vec, const VectorType<_T, n>& other)
    {
        VectorType<_T, n> res;
        for(uint16_t i = 0; i < n; i++)
            res.elements[i] = vec.elements[i] + other.elements[i];
        return res;
    }

    template<typename _T, uint16_t n>
    static VectorType<_T, n> Subtract(const VectorType<_T, n>& vec, const VectorType<_T, n>& other)
    {
        VectorType<_T, n> res;
        for(uint16_t i = 0; i < n; i++)
            res.elements[i] = vec.elements[i] - other.elements[i];
        return res;
    }

    template<typename _T, uint16_t n>
    static VectorType<_T, n> Multiply(const VectorType<_T, n>& vec, const float scalar)
    {
        VectorType<_T, n> res;
        for(uint16_t i = 0; i < n; i++)
            res.elements[i] = vec.elements[i] * scalar;
        return res;
    }

    template<typename _T, uint16_t n>
    static VectorType<_T, n> Divide(const VectorType<_T, n>& vec, const float scalar)
    {
        VectorType<_T, n> res;
        for(uint16_t i = 0; i < n; i++)
            res.elements[i] = vec.elements[i] / scalar;
        return res;
    }

    template<typename _T, uint16_t n>
    static VectorType<_T, n> operator+(const VectorType<_T, n>& left, const VectorType<_T, n>& right)
    {
        return Add(left, right);
    }

    template<typename _T, uint16_t n>
    static VectorType<_T, n> operator-(const VectorType<_T, n>& left, const VectorType<_T, n>& right)
    {
        return Subtract(left, right);
    }

    template<typename _T, uint16_t n>
    static VectorType<_T, n> operator+=(VectorType<_T, n>& left, const VectorType<_T, n>& right)
    {
        for(uint16_t i = 0; i < n; i++)
            left.elements[i] += right.elements[i];
        return left;
    }

    template<typename _T, uint16_t n>
    static VectorType<_T, n> operator-=(VectorType<_T, n>& left, const VectorType<_T, n>& right)
    {
        for(uint16_t i = 0; i < n; i++)
            left.elements[i] -= right.elements[i];
        return left;
    }

    template<typename _T, uint16_t n>
    static VectorType<_T, n> operator*(const VectorType<_T, n>& left, const float right)
    {
        return Multiply(left, right);
    }

    template<typename _T, uint16_t n>
    static VectorType<_T, n> operator/(const VectorType<_T, n>& left, const float right)
    {
        return Divide(left, right);
    }

    template<typename _T, uint16_t n>
    static VectorType<_T, n> operator*=(VectorType<_T, n>& left, const float right)
    {
        for(uint16_t i = 0; i < n; i++)
            left.elements[i] *= right;
        return left;
    }

    template<typename _T, uint16_t n>
    static VectorType<_T, n> operator/=(VectorType<_T, n>& left, const float right)
    {
        for(uint16_t i = 0; i < n; i++)
            left.elements[i] /= right;
        return left;
    }

    template<typename _T, uint16_t n>
    static bool operator==(const VectorType<_T, n>& left, const VectorType<_T, n>& right)
    {
        for(uint16_t i = 0; i < n; i++)
        {
            if(left.elements[i] != right.elements[i])
                return false;
        }
        return true;
    }

    template<typename _T, uint16_t n>
    static bool operator!=(const VectorType<_T, n>& left, const VectorType<_T, n>& right)
    {
        return !(left == right);
    }

    template<typename _T, uint16_t n>
    static VectorType<_T, n> operator-(VectorType<_T, n> vec)
    {
        return vec.Negate();
    }

    template<typename _T, uint16_t n>
    static std::ostream& operator<<(std::ostream& stream, const VectorType<_T, n>& vec)
    {
        stream << "( ";
        for(int i = 0; i < 4; i++)
        {
            stream<< vec.elements[i];
            if(i != 3)
                stream<< ", ";
        }
        stream<< ")\n";
        return stream;
    }

    template<typename _T, uint16_t n>
    VectorType<_T, n> Normalize(const VectorType<_T, n>& vec)
    {
        return vec / vec.Magnitude();
    }

    template<typename _T, uint16_t n>
    GMATH_TYPE_PRECISION Dot(const VectorType<_T, n>& vec, const VectorType<_T, n>& other)
    {
        GMATH_TYPE_PRECISION total = 0.0f;
        for(uint16_t i = 0; i < n; i++)
            total += vec.elements[i] * other.elements[i];
        return sqrt(total);
    }

    template <typename _T>
    VectorType<_T, 3> Cross(const VectorType<_T, 3>& vec, const VectorType<_T, 3>& other)
    {
        VectorType<_T, 3> res;
        res.elements[0] = vec.elements[1] * other.elements[2] - other.elements[1] * vec.elements[2];
        res.elements[1] = vec.elements[0] * other.elements[2] - other.elements[0] * vec.elements[2];
        res.elements[2] = vec.elements[0] * other.elements[1] - other.elements[0] * vec.elements[1];
        return res;
    }

    template<uint16_t n>
    using Vector = VectorType<GMATH_TYPE_PRECISION, n>;

    using Vector4 = Vector<4>;
    using Vector3 = Vector<3>;
    using Vector2 = Vector<2>;

}

#endif // VECTOR_GMATH_H
