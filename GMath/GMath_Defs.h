#ifndef LOGS_GMATH_H
#define LOGS_GMATH_H

#define uint16_t unsigned short
#define uint32_t unsigned int

#define M_PI		3.14159265358979323846

namespace GMath{

    static inline float toRadians(float degrees)
    {
        return degrees * (M_PI / 180.0f);
    }

    static inline float toDegrees(float radians)
    {
        return radians * (180.0f/ M_PI);
    }

}

#define GMATH_TYPE_PRECISION float

#endif // LOGS_GMATH_H
