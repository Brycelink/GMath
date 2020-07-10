#ifndef QUATERNION_AXIS_H
#define QUATERNION_AXIS_H

#include "Logs.h"
#include "GMath/Vector.h"

namespace GMath{

    template<typename _T>
    struct QuaternionType
    {
    public:
        union{
            struct{_T x, y, z, w;};
            _T elements[4];
        };

    public:
        QuaternionType()
        {
            x = 0.0f;
            y = 0.0f;
            z = 0.0f;
            w = 0.0f;
        }

        QuaternionType(const QuaternionType<_T>& other)
        {
            x = other.x;
            y = other.y;
            z = other.z;
            w = other.w;
        }

        QuaternionType(const _T& _x, const _T& _y, const _T& _z, const _T& _w)
        {
            x = _x;
            y = _y;
            z = _z;
            w = _w;
        }

        QuaternionType(_T data[4])
        {
            x = data[0];
            y = data[1];
            z = data[2];
            w = data[3];
        }

        GMATH_TYPE_PRECISION Magnitude()
        {
            return sqrt((x * x) + (y * y) + (z * z) + (w * w));
        }

        QuaternionType<_T>& Normalize()
        {
            GMATH_TYPE_PRECISION length = (*this).Magnitude();

            x /= length;
            y /= length;
            z /= length;
            w /= length;

            return (*this);
        }

        QuaternionType<_T> Conjugate()
        {
            return Quaternion(-x, -y, -z, w);
        }

        QuaternionType<_T>& Multiply(const QuaternionType<_T>& other)
        {
            float _w = w * other.w - x * other.x - y * other.y - z * other.z;
            float _x = x * other.w + w * other.x + y * other.z - z * other.y;
            float _y = y * other.w + w * other.y + z * other.x - x * other.z;
            float _z = z * other.w + w * other.z + x * other.y - y * other.x;

            x = _x;
            y = _y;
            z = _z;
            w = _w;

            return (*this);
        }

        QuaternionType<_T>& Multiply(const Vector3& vec)
        {
            float _w = -x * vec.x - y * vec.y - z * vec.z;
            float _x =  w * vec.x + y * vec.z - z * vec.y;
            float _y =  w * vec.y + z * vec.x - x * vec.z;
            float _z =  w * vec.z + x * vec.y - y * vec.x;

            x = _x;
            y = _y;
            z = _z;
            w = _w;

            return (*this);
        }

        friend QuaternionType<_T> operator*(QuaternionType left, const Vector3& right)
        {
            return left.Multiply(right);
        }

        friend QuaternionType operator*(QuaternionType<_T> left, const QuaternionType<_T>& right)
        {
            return left.Multiply(right);
        }

        bool operator==(const QuaternionType<_T>& other)
        {
            return ((x == other.x) && (y == other.y) && (z == other.z) && (w == other.w));
        }

        bool operator!=(const QuaternionType<_T>& other)
        {
            return !((*this) == other);
        }

        QuaternionType<_T>& operator*=(const QuaternionType<_T>& other)
        {
            return ((*this).Multiply(other));
        }

        friend std::ostream& operator<<(std::ostream &stream, const QuaternionType<_T>& quart)
        {
            stream << "( " << quart.x << ", " << quart.y << ", " << quart.z << ", " << quart.w << ")";
            return stream;
        }

        void Print()
        {
            printf("(%f, %f, %f, %f)", x, y, z, w);
        }

    };

    using Quaternion = QuaternionType<GMATH_TYPE_PRECISION>;

}

#endif // QUATERNION_AXIS_H
