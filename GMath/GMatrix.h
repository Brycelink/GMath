#ifndef G_MATRIX_GMATH_H
#define G_MATRIX_GMATH_H

#include "GMath/Matrix.h"

namespace GMath{

    template<uint16_t rows, uint16_t cols>
    using MatrixC = MatrixTypeC<GMATH_TYPE_PRECISION, rows, cols>;

    template<uint16_t n>
    using MatrixSqC = MatrixC<n, n>;

    template<uint16_t rows, uint16_t cols>
    using MatrixR = MatrixTypeR<GMATH_TYPE_PRECISION, rows, cols>;

    template<uint16_t n>
    using MatrixSqR = MatrixR<n, n>;

    template<uint16_t rows, uint16_t cols>
    using Matrix = MatrixC<rows, cols>;

    template<uint16_t n>
    using MatrixSq = Matrix< n, n>;

    typedef MatrixSq<4> Matrix4;
    typedef MatrixSq<3> Matrix3;
    typedef MatrixSq<2> Matrix2;
    typedef MatrixSq<1> Matrix1;

}

#endif // G_MATRIX_GMATH_H
