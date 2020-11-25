#ifndef GMATH_GMATH_H
#define GMATH_GMATH_H

#include "GMath/GMatrix.h"
#include "GMath/GVector.h"
#include "GMath/GQuaternion.h"

//GMath supports use of both column Major and Row Major for Matrices
/*
N/B:    (i)   MatrixTypeC represents column major matrices
        (ii)  MatrixTypeR represents row major matrices
        (iii) GMATH_TYPE_PRECISION is set to float by default

The Code below demonstrates simplifying GMath classes to more meaningful
and easier to handle syntax as specified in the GMatrix, GVector and GQuaternion
header files. It specifies use of Column Major for matrices by default
and use of floating point values for the elements of vectors and matrices.
Modification of GMatrix, GVector and GQuartenion to suite user's needs is
encouraged.

---------------------------------------------------------------------------

namespace GMath
{
    //----------GMatrix.h-----------------------------------------
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

    //----------GVector.h-----------------------------------------
    template<uint16_t n>
    using Vector = VectorType<GMATH_TYPE_PRECISION, n>;

    using Vector4 = Vector<4>;
    using Vector3 = Vector<3>;
    using Vector2 = Vector<2>;

    //----------GQuaternion.h-----------------------------------------

    using Quaternion = QuaternionType<GMATH_TYPE_PRECISION>;

}

--------------------------------------------------------------------------

*/

#endif // GMATH_GMATH_H
