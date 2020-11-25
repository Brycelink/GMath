#ifndef G_VECTOR_GMATH_H
#define G_VECTOR_GMATH_H

#include "GMath/Vector.h"

namespace GMath{

    template<uint16_t n>
    using Vector = VectorType<GMATH_TYPE_PRECISION, n>;

    using Vector4 = Vector<4>;
    using Vector3 = Vector<3>;
    using Vector2 = Vector<2>;

}

#endif // G_VECTOR_GMATH_H
